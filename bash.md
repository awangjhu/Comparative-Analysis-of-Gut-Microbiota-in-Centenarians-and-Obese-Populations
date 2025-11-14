# BUDDY Project: Microbiome Analysis with Kraken2

# Commands

## Setup Instructions
conda activate qb25

## Download paired-end reads
wget -i qbshortened_download_data.txt --continue #qbshortened_download_data.txt contains
### qbshortened_download_data.txt
This text file contains the FTP download links for the paired-end FASTQ files used in this project.  
Each line corresponds to a sample (e.g., centenarian `C_7G6` or young `Y_7G39`), listing the URLs for the forward and reverse reads to be fetched using the `wget` command.

## Install latest version of kraken software
conda install -c bioconda kraken2=2.1.6

## Create new directory for the Kraken Database
in /BUDDYproject --> mkdir krakendb

## Install minikraken 8GB database (good for working at home since the actual Kraken database is 100GB) done in the krakendb subdirectory
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz

## Extract MiniKraken2 database
tar -xvzf minikraken2_v2_8GB_201904.tgz

## Decompress Raw Data files (ran from cd BUDDYproject)
gzip -dc S_7G6_1.fq.gz  > rawreads/S_7G6_1.fq
gzip -dc S_7G6_2.fq.gz  > rawreads/S_7G6_2.fq
gzip -dc S_7G39_1.fq.gz > rawreads/S_7G39_1.fq
gzip -dc S_7G39_2.fq.gz > rawreads/S_7G39_2.fq

## Running Kraken2 Classification
# The commands below classify each read pair against the MiniKraken2 database using 6 threads for faster processing while keeping the system responsive. The --report flag creates a summary file with taxonomic abundances, while > results/..._output.txt stores individual read classifications.
# Centenarian sample
kraken2 \
  --db krakendb/minikraken2_v2_8GB_201904_UPDATE \
  --paired \
  --threads 6 \
  --use-names \
  --report results/centenarian_S_7G6_report.txt \
  rawreads/S_7G6_1.fq rawreads/S_7G6_2.fq \
  > results/centenarian_S_7G6_output.txt


# Young sample 
kraken2 \
  --db krakendb/minikraken2_v2_8GB_201904_UPDATE \
  --paired \
  --threads 6 \
  --use-names \
  --report results/young_S_7G39_report.txt \
  rawreads/S_7G39_1.fq rawreads/S_7G39_2.fq \
  > results/young_S_7G39_output.txt

  # Convert the report.txt files to .csv for R analysis
awk -F'\t' 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6}' centenarian_S_7G6_report.txt > centenarian_S_7G6_report.csv
awk -F'\t' 'BEGIN{OFS=","} {print $1,$2,$3,$4,$5,$6}' young_S_7G39_report.txt > young_S_7G39_report.csv