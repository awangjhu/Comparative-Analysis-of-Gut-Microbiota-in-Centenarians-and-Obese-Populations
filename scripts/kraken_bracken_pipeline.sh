#!/usr/bin/env bash
# kraken_bracken_pipeline.sh - Fixed for 2020 Kraken2 database (Newer data base had issues :( ))
# In summary: this runs Kraken2 → Fixes rank codes → Bracken on paired FASTQs

set -euo pipefail
IFS=$'\n\t'

DB=""
READS_DIR="./rawreads"
OUT_DIR="./results"
THREADS=8
READ_LEN=150
LEVELS="S,G,P"
FORCE=0
KRAKEN_OPTS=""
BRACKEN_THRESH=10

log(){ echo -e "[\e[1;34mINFO\e[0m] $*"; }
warn(){ echo -e "[\e[1;33mWARN\e[0m] $*" >&2; }
err(){ echo -e "[\e[1;31mERR \e[0m] $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || err "Required tool '$1' not found in PATH."; }
abspath(){ perl -MCwd -e 'print Cwd::abs_path(shift)' "$1"; }

print_help(){
  cat <<'H'
Usage:
  bash kraken_bracken_pipeline.sh \
    --db /path/to/kraken_db \
    --reads ./rawreads \
    --out ./results \
    --threads 8 \
    --levels S,G,P \
    --read-len 150
H
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) print_help; exit 0 ;;
    --db) DB=$(abspath "$2"); shift 2 ;;
    --reads) READS_DIR=$(abspath "$2"); shift 2 ;;
    --out) OUT_DIR=$(abspath "$2"); shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --read-len) READ_LEN="$2"; shift 2 ;;
    --levels) LEVELS="$2"; shift 2 ;;
    --force) FORCE="$2"; shift 2 ;;
    --kraken-opts) KRAKEN_OPTS="$2"; shift 2 ;;
    --bracken-thresh) BRACKEN_THRESH="$2"; shift 2 ;;
    *) err "Unknown arg: $1" ;;
  esac
done

[[ -n "$DB" ]] || err "--db is required"
[[ -d "$READS_DIR" ]] || err "--reads directory not found: $READS_DIR"

need kraken2
need bracken

# Use Bracken v2.9 
if [[ -x "$HOME/Bracken_v2.9/bracken" ]]; then
  BRACKEN_CMD="$HOME/Bracken_v2.9/bracken"
  log "Using Bracken v2.9 from $HOME/Bracken_v2.9"
else
  BRACKEN_CMD="bracken"
  warn "Using system bracken - v2.9 recommended for this database"
fi

command -v pigz >/dev/null 2>&1 || warn "pigz not found; gzip may be slower"
command -v parallel >/dev/null 2>&1 || warn "GNU parallel not found; running serially"

# DIR
REPORT_DIR="$OUT_DIR/reports"
mkdir -p "$REPORT_DIR"
IFS=',' read -r -a LEVEL_ARR <<< "$LEVELS"
for L in "${LEVEL_ARR[@]}"; do
  mkdir -p "$OUT_DIR/bracken/$L"
done

# Find Pairs
find_pairs(){
  shopt -s nullglob
  local exts=("fq.gz" "fastq.gz" "fq" "fastq")
  local pairs=()
  local count=0
  for ext in "${exts[@]}"; do
    for R1 in "$READS_DIR"/*_1.$ext; do
      local R2="${R1/_1.$ext/_2.$ext}"
      [[ -f "$R2" ]] || { warn "Missing mate for $R1"; continue; }
      local sample
      sample=$(basename "$R1")
      sample="${sample%_1.$ext}"
      pairs+=("$(printf '%s\t%s\t%s' "$sample" "$R1" "$R2")")
      count=$((count+1))
    done
  done
  printf '%s\n' "${pairs[@]}"
}

SAMPLES_TSV=$(mktemp)
find_pairs > "$SAMPLES_TSV"
[[ -s "$SAMPLES_TSV" ]] || err "No paired FASTQ files found in $READS_DIR"
log "Discovered $(wc -l < "$SAMPLES_TSV") paired samples."

# Kraken2
run_kraken2(){
  local sample="$1" R1="$2" R2="$3"
  local krep="$REPORT_DIR/${sample}.kreport"
  if [[ -f "$krep" && "$FORCE" -ne 1 ]]; then
    log "(skip) Kraken2 report exists: $krep"
    return 0
  fi
  log "Kraken2 → $sample"
  kraken2_args=(
    --db "$DB"
    --paired "$R1" "$R2"
    --threads "$THREADS"
    --report "${krep}"
    --output /dev/null
  )
  if [[ -n "${KRAKEN_OPTS}" ]]; then
    extra=( ${KRAKEN_OPTS} )
    kraken2_args+=( "${extra[@]}" )
  fi
  kraken2 "${kraken2_args[@]}"
}

export -f run_kraken2
export REPORT_DIR DB THREADS FORCE KRAKEN_OPTS

while IFS=$'\t' read -r sample R1 R2; do
  run_kraken2 "$sample" "$R1" "$R2"
done < "$SAMPLES_TSV"

# Braken
log "Fixing Kraken2 reports for Bracken compatibility..."
for krep in "$REPORT_DIR"/*.kreport; do
  [[ -f "$krep" ]] || continue
  python3 - "$krep" "$krep.tmp" << 'PYEOF'
import sys
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
for i, line in enumerate(lines):
    parts = line.split('\t')
    if len(parts) >= 5:
        rank = parts[3].strip()
        if rank.isdigit() or rank == '':
            taxid = parts[4].strip()
            if taxid in ['131567', '2', '2157', '2759', '10239', '10442']:
                parts[3] = 'D'
            else:
                parts[3] = 'U1'
            lines[i] = '\t'.join(parts)
with open(sys.argv[2], 'w') as f:
    f.writelines(lines)
PYEOF
  mv "$krep.tmp" "$krep"
done
log "Reports fixed successfully"

run_bracken(){
  local sample="$1" level="$2"
  local krep="$REPORT_DIR/${sample}.kreport"
  local outdir="$OUT_DIR/bracken/$level"
  local bout="$outdir/${sample}.bracken"
  local brep="$outdir/${sample}.bracken.report.txt"

  [[ -s "$krep" ]] || err "Missing Kraken2 report for $sample: $krep"
  mkdir -p "$outdir"

  if [[ -f "$bout" && "$FORCE" -ne 1 ]]; then
    log "(skip) Bracken ($level) exists: $bout"
    return 0
  fi

  # Convert to RELATIVE paths for Bracken
  local krep_rel bout_rel brep_rel
  krep_rel="$(python -c 'import os,sys; print(os.path.relpath(sys.argv[1], os.getcwd()))' "$krep")"
  bout_rel="$(python -c 'import os,sys; print(os.path.relpath(sys.argv[1], os.getcwd()))' "$bout")"
  brep_rel="$(python -c 'import os,sys; print(os.path.relpath(sys.argv[1], os.getcwd()))' "$brep")"

  log "Bracken ($level) → $sample"
  $BRACKEN_CMD \
    -d "$DB" \
    -i "$krep_rel" \
    -o "$bout_rel" \
    -w "$brep_rel" \
    -r "$READ_LEN" \
    -l "$level" \
    -t "$BRACKEN_THRESH"
}

export -f run_bracken
export OUT_DIR READ_LEN BRACKEN_THRESH BRACKEN_CMD

while IFS=$'\t' read -r sample _R1 _R2; do
  for L in "${LEVEL_ARR[@]}"; do
    run_bracken "$sample" "$L"
  done
done < "$SAMPLES_TSV"

rm -f "$SAMPLES_TSV"

log "Done. See:"
for L in "${LEVEL_ARR[@]}"; do
  log "  - $OUT_DIR/bracken/$L/<SAMPLE>.bracken"
  log "  - $OUT_DIR/bracken/$L/<SAMPLE>.bracken.report.txt"
done