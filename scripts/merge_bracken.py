#!/usr/bin/env python3
# Merge Bracken results for Sardinia centenarian study.
# Designed to work with kraken_bracken_pipeline.sh output structure.
# BUDDY MERGER!!!

import os
import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict


def load_metadata(metadata_file):
    # Load sample metadata file.
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')
        required_cols = ['sample_id', 'age_group']
        for col in required_cols:
            if col not in metadata.columns:
                print(f"Error: metadata file must have '{col}' column", file=sys.stderr)
                sys.exit(1)
        print(f" Loaded metadata for {len(metadata)} samples")
        print(f"  Groups: {', '.join(metadata['age_group'].unique())}")
        return metadata
    except Exception as e:
        print(f"Metadata loading issue, need to fix: {e}", file=sys.stderr)
        sys.exit(1)


def parse_bracken_file(filepath):
    #Parse a single Bracken output file
    try:
        # Bracken output format:
        df = pd.read_csv(filepath, sep='\t')
        if 'name' not in df.columns or 'new_est_reads' not in df.columns:
            print(f"Expected columns not found in {filepath}", file=sys.stderr)
            return None
        return df
    except Exception as e:
        print(f"Could not parse {filepath}: {e}", file=sys.stderr)
        return None


def merge_by_level(bracken_dir, output_dir, levels, metadata=None):
    
    # Merge all Bracken files for each taxonomic level into wide-format tables.
    # Creates abundance tables, relative abundance tables, and group averages.
    
    bracken_path = Path(bracken_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    level_names = {
        'S': 'species',
        'G': 'genus', 
        'P': 'phylum',
        'C': 'class',
        'O': 'order',
        'F': 'family'
    }
    
    # Create sample to group mapping
    sample_groups = {}
    if metadata is not None:
        sample_groups = dict(zip(metadata['sample_id'], metadata['age_group']))
    
    print("\n" + "="*70)
    print("MERGING BRACKEN RESULTS")
    print("="*70)
    
    for level in levels:
        level_dir = bracken_path / level
        if not level_dir.exists():
            print(f"\n Directory not found: {level_dir}")
            continue
        
        level_name = level_names.get(level, level)
        print(f"\n Processing taxonomic level: {level} ({level_name.upper()})")
        print("-" * 70)
        
        # Collect all samples for this level
        all_data = {}
        sample_metadata_map = {}
        
        bracken_files = list(level_dir.glob("*.bracken"))
        print(f"   Found {len(bracken_files)} Bracken files")
        
        for bracken_file in bracken_files:
            sample_name = bracken_file.stem  # S_7G38 (removes .bracken)
            df = parse_bracken_file(bracken_file)
            
            if df is not None and not df.empty:
                # Extract taxonomy name and abundance
                taxa_abundance = dict(zip(df['name'], df['new_est_reads']))
                all_data[sample_name] = taxa_abundance
                
                # Store metadata for this sample
                if sample_name in sample_groups:
                    sample_metadata_map[sample_name] = sample_groups[sample_name]
                else:
                    sample_metadata_map[sample_name] = 'unknown'
        
        if not all_data:
            print(f" No data found for level {level}")
            continue
        
        print(f" Successfully parsed {len(all_data)} samples")
        
        # Convert to format DataFrame (rows=taxa, columns=samples)
        abundance_df = pd.DataFrame(all_data).fillna(0)
        
        # Sort columns by group if metadata available
        if sample_metadata_map:
            groups_order = ['young', 'elderly', 'centenarian', 'unknown']
            col_order = sorted(
                abundance_df.columns, 
                key=lambda x: (
                    groups_order.index(sample_metadata_map.get(x, 'unknown')) 
                    if sample_metadata_map.get(x, 'unknown') in groups_order 
                    else 999,
                    x
                )
            )
            abundance_df = abundance_df[col_order]
        
        # Sort rows by total abundance (most abundant taxa first)
        abundance_df['_total'] = abundance_df.sum(axis=1)
        abundance_df = abundance_df.sort_values('_total', ascending=False)
        abundance_df = abundance_df.drop('_total', axis=1)
        
        print(f"   Taxa detected: {len(abundance_df)}")
        print(f"   Most abundant: {abundance_df.sum(axis=1).idxmax()}")
        
        # Output 1: Raw abundance table
        output_file = output_path / f"{level_name}_abundance.tsv"
        abundance_df.to_csv(output_file, sep='\t')
        print(f"\n {output_file.name}")
        print(f" ({len(abundance_df)} taxa × {len(abundance_df.columns)} samples)")
        
        # Output 2: Relative abundance
        rel_abundance_df = abundance_df.div(abundance_df.sum(axis=0), axis=1) * 100
        rel_output_file = output_path / f"{level_name}_relative_abundance.tsv"
        rel_abundance_df.to_csv(rel_output_file, sep='\t')
        print(f" {rel_output_file.name}")
        print(f" (relative abundance %)")
        
        # Output 3: Annotations
        if sample_metadata_map and any(v != 'unknown' for v in sample_metadata_map.values()):
            annotated_output = output_path / f"{level_name}_abundance_annotated.tsv"
            with open(annotated_output, 'w') as f:
                # Write header with sample names
                f.write('taxon\t' + '\t'.join(abundance_df.columns) + '\n')
                # Write group labels
                f.write('GROUP\t' + '\t'.join([
                    sample_metadata_map.get(col, 'unknown') 
                    for col in abundance_df.columns
                ]) + '\n')
                # Write data
                abundance_df.to_csv(f, sep='\t', header=False)
            print(f" {annotated_output.name}")
            print(f"(with group labels in header)")
        
        # Output 4: Group analysis
        if sample_metadata_map and any(v != 'unknown' for v in sample_metadata_map.values()):
            # Organize samples by group
            groups = defaultdict(list)
            for sample, group in sample_metadata_map.items():
                if group != 'unknown' and sample in abundance_df.columns:
                    groups[group].append(sample)
            
            if len(groups) > 0:
                # Calculate mean and std for each group
                group_stats = pd.DataFrame()
                
                # Sort groups in logical order
                group_order = ['young', 'elderly', 'centenarian']
                for group_name in group_order:
                    if group_name in groups:
                        samples = groups[group_name]
                        group_stats[f'{group_name}_mean'] = abundance_df[samples].mean(axis=1)
                        group_stats[f'{group_name}_std'] = abundance_df[samples].std(axis=1)
                        group_stats[f'{group_name}_n'] = len(samples)
                
                # Add relative abundance version
                group_rel_stats = pd.DataFrame()
                for group_name in group_order:
                    if group_name in groups:
                        samples = groups[group_name]
                        group_rel_stats[f'{group_name}_mean'] = rel_abundance_df[samples].mean(axis=1)
                        group_rel_stats[f'{group_name}_std'] = rel_abundance_df[samples].std(axis=1)
                
                # Save absolute abundance group averages
                group_output = output_path / f"{level_name}_group_averages.tsv"
                group_stats.to_csv(group_output, sep='\t')
                print(f"{group_output.name}")
                print(f"      (group means & std devs ")
                
                # Save relative abundance group averages
                group_rel_output = output_path / f"{level_name}_group_relative_averages.tsv"
                group_rel_stats.to_csv(group_rel_output, sep='\t')
                print(f" {group_rel_output.name}")
                print(f" (group means in %")
                
                # Print group sample counts
                print(f"\n Sample counts per group:")
                for group_name in group_order:
                    if group_name in groups:
                        print(f" {group_name}: n={len(groups[group_name])}")


def create_sample_summary(bracken_dir, output_dir, levels, metadata=None):
    
    # Summary table showing read counts and taxa detected per sample.
    
    bracken_path = Path(bracken_dir)
    output_path = Path(output_dir)
    
    sample_info = defaultdict(dict)
    
    for level in levels:
        level_dir = bracken_path / level
        if not level_dir.exists():
            continue
            
        for bracken_file in level_dir.glob("*.bracken"):
            sample_name = bracken_file.stem
            df = parse_bracken_file(bracken_file)
            
            if df is not None and 'new_est_reads' in df.columns:
                total_reads = df['new_est_reads'].sum()
                num_taxa = len(df)
                sample_info[sample_name][f'{level}_reads'] = int(total_reads)
                sample_info[sample_name][f'{level}_taxa'] = num_taxa
    
    if sample_info:
        summary_df = pd.DataFrame.from_dict(sample_info, orient='index')
        summary_df.index.name = 'sample_id'
        
        # Add metadata
        if metadata is not None:
            metadata_indexed = metadata.set_index('sample_id')
            summary_df = summary_df.join(metadata_indexed, how='left')
            if 'age_group' in summary_df.columns:
                cols = ['age_group'] + [col for col in summary_df.columns if col != 'age_group']
                summary_df = summary_df[cols]
        
        summary_file = output_path / "sample_summary.tsv"
        summary_df.to_csv(summary_file, sep='\t')
        
        print("\n" + "="*70)
        print(f" Sample summary saved: {summary_file}")
        print("="*70)
        
        # Print some statistics
        if 'age_group' in summary_df.columns:
            print("\nSamples per group:")
            print(summary_df['age_group'].value_counts().to_string())


def print_analysis_guide(output_dir):
    """Print helpful guidance for next steps."""
    print("\n" + "="*70)
    print("BUDDY MERGE COMPLETE!")
    print("="*70)
    print("\n Output files in:", output_dir)
    print("\n KEY FILES FOR ANALYSIS:")
    print("   1. species_group_averages.tsv")
    print("     comparing centenarian vs young vs elderly")
    print("   2. species_relative_abundance.tsv")
    print("      individual sample analysis")
    print("   3. sample_summary.tsv")
    print("      QC: check read counts and taxa detected")

def main():
    parser = argparse.ArgumentParser(
        description='Merge Bracken results for Sardinia centenarian microbiome study',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        """
    )
    parser.add_argument('--input', '-i', required=True,
                       help='Input directory containing Bracken results (with S/, G/, P/ subdirs)')
    parser.add_argument('--output', '-o', required=True,
                       help='Output directory for merged tables')
    parser.add_argument('--levels', '-l', default='S,G,P',
                       help='Comma-separated taxonomic levels')
    parser.add_argument('--metadata', '-m', 
                       help='Metadata file with sample_id and age_group columns')
    
    args = parser.parse_args()
    
    # Parse levels
    levels = [l.strip() for l in args.levels.split(',')]
    
    # Load metadata 
    metadata = None
    if args.metadata:
        if not os.path.exists(args.metadata):
            print(f"Metadata file not found: {args.metadata}", file=sys.stderr)
            sys.exit(1)
        metadata = load_metadata(args.metadata)
    else:
        print("NO METADATA BUDDY")
        print(" ADD: metadata sardinia_sample_metadata.txt")
    
    print(f"\n Input directory: {args.input}")
    print(f" Output directory: {args.output}")
    print(f" Processing levels: {', '.join(levels)}")
    
    # Input directory
    if not os.path.exists(args.input):
        print(f"\n Input directory is not found, buddy: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Merge
    merge_by_level(args.input, args.output, levels, metadata)
    create_sample_summary(args.input, args.output, levels, metadata)

    print_analysis_guide(args.output)