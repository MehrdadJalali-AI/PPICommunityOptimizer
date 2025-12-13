#!/usr/bin/env python3
"""
MCL vs LEAF-PPI Comparison Analysis

This script creates a focused comparison between MCL (baseline) and LEAF-PPI (optimized)
to demonstrate the improvement in community detection quality.
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
import json

def load_comparison_data():
    """Load comparison data from CSV files."""
    try:
        df_string = pd.read_csv('results_string_updated.csv')
        df_gavin = pd.read_csv('results_gavin_updated.csv')
        return df_string, df_gavin
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please run 'python generate_updated_results.py' first to generate comparison data.")
        sys.exit(1)

def extract_mcl_vs_leaf(df, dataset_name):
    """Extract MCL and LEAF-PPI results from dataframe."""
    mcl_row = df[df['method'] == 'MCL'].iloc[0] if len(df[df['method'] == 'MCL']) > 0 else None
    leaf_row = df[df['method'] == 'LEA_Overlapping'].iloc[0] if len(df[df['method'] == 'LEA_Overlapping']) > 0 else None
    
    if mcl_row is None or leaf_row is None:
        print(f"Warning: Missing MCL or LEAF-PPI data for {dataset_name}")
        return None
    
    return {
        'dataset': dataset_name,
        'mcl': mcl_row.to_dict(),
        'leaf': leaf_row.to_dict()
    }

def calculate_improvement(mcl_value, leaf_value, metric_name, higher_is_better=True):
    """Calculate improvement percentage."""
    if pd.isna(mcl_value) or pd.isna(leaf_value):
        return None, None
    
    if mcl_value == 0:
        if leaf_value == 0:
            return 0.0, "No change"
        return float('inf'), "∞% improvement (from 0)"
    
    if higher_is_better:
        improvement = ((leaf_value - mcl_value) / abs(mcl_value)) * 100
    else:
        improvement = ((mcl_value - leaf_value) / abs(mcl_value)) * 100
    
    if improvement > 0:
        direction = "improvement"
    elif improvement < 0:
        direction = "degradation"
    else:
        direction = "no change"
    
    return improvement, direction

def create_comparison_table(data):
    """Create a detailed comparison table."""
    metrics = [
        ('modularity', 'Modularity', True, 'Higher is better'),
        ('conductance', 'Conductance', False, 'Lower is better'),
        ('num_communities', 'Number of Communities', None, 'Context-dependent'),
        ('avg_community_size', 'Avg Community Size', None, 'Context-dependent'),
        ('mean_go_jaccard', 'GO Jaccard Similarity', True, 'Higher is better'),
        ('runtime_sec', 'Runtime (seconds)', False, 'Lower is better'),
    ]
    
    rows = []
    for metric_key, metric_name, higher_is_better, note in metrics:
        row = {
            'Metric': metric_name,
            'Note': note,
        }
        
        for dataset_data in data:
            if dataset_data is None:
                continue
            
            dataset_name = dataset_data['dataset']
            mcl_val = dataset_data['mcl'].get(metric_key)
            leaf_val = dataset_data['leaf'].get(metric_key)
            
            row[f'{dataset_name} - MCL'] = f"{mcl_val:.6f}" if pd.notna(mcl_val) else "N/A"
            row[f'{dataset_name} - LEAF-PPI'] = f"{leaf_val:.6f}" if pd.notna(leaf_val) else "N/A"
            
            if higher_is_better is not None:
                improvement, direction = calculate_improvement(mcl_val, leaf_val, metric_name, higher_is_better)
                if improvement is not None:
                    row[f'{dataset_name} - Improvement'] = f"{improvement:+.2f}% ({direction})"
                else:
                    row[f'{dataset_name} - Improvement'] = "N/A"
        
        rows.append(row)
    
    return pd.DataFrame(rows)

def create_summary_statistics(data):
    """Create summary statistics showing overall improvement."""
    summary = {
        'dataset': [],
        'metric': [],
        'mcl_value': [],
        'leaf_value': [],
        'improvement_pct': [],
        'improvement_direction': []
    }
    
    key_metrics = [
        ('modularity', True),
        ('conductance', False),
        ('mean_go_jaccard', True),
    ]
    
    for dataset_data in data:
        if dataset_data is None:
            continue
        
        dataset_name = dataset_data['dataset']
        
        for metric_key, higher_is_better in key_metrics:
            mcl_val = dataset_data['mcl'].get(metric_key)
            leaf_val = dataset_data['leaf'].get(metric_key)
            
            if pd.isna(mcl_val) or pd.isna(leaf_val):
                continue
            
            improvement, direction = calculate_improvement(mcl_val, leaf_val, metric_key, higher_is_better)
            
            summary['dataset'].append(dataset_name)
            summary['metric'].append(metric_key)
            summary['mcl_value'].append(mcl_val)
            summary['leaf_value'].append(leaf_val)
            summary['improvement_pct'].append(improvement if improvement is not None else 0)
            summary['improvement_direction'].append(direction)
    
    return pd.DataFrame(summary)

def print_comparison_report(data, comparison_df, summary_df):
    """Print a formatted comparison report."""
    print("=" * 100)
    print("MCL vs LEAF-PPI: Community Detection Quality Comparison")
    print("=" * 100)
    print()
    
    print("SUMMARY OF IMPROVEMENTS")
    print("-" * 100)
    print(summary_df.to_string(index=False))
    print()
    
    print("DETAILED METRIC COMPARISON")
    print("-" * 100)
    print(comparison_df.to_string(index=False))
    print()
    
    # Key findings
    print("KEY FINDINGS")
    print("-" * 100)
    
    for dataset_data in data:
        if dataset_data is None:
            continue
        
        dataset_name = dataset_data['dataset']
        mcl = dataset_data['mcl']
        leaf = dataset_data['leaf']
        
        print(f"\n{dataset_name} Dataset:")
        print(f"  • Modularity: {mcl.get('modularity', 'N/A'):.6f} → {leaf.get('modularity', 'N/A'):.6f}")
        if pd.notna(mcl.get('modularity')) and pd.notna(leaf.get('modularity')):
            mod_improvement, _ = calculate_improvement(mcl.get('modularity'), leaf.get('modularity'), 'modularity', True)
            if mod_improvement is not None:
                print(f"    Improvement: {mod_improvement:+.2f}%")
        
        print(f"  • Conductance: {mcl.get('conductance', 'N/A'):.6f} → {leaf.get('conductance', 'N/A'):.6f}")
        if pd.notna(mcl.get('conductance')) and pd.notna(leaf.get('conductance')):
            cond_improvement, _ = calculate_improvement(mcl.get('conductance'), leaf.get('conductance'), 'conductance', False)
            if cond_improvement is not None:
                print(f"    Improvement: {cond_improvement:+.2f}%")
        
        print(f"  • Communities: {mcl.get('num_communities', 'N/A')} → {leaf.get('num_communities', 'N/A')}")
        print(f"  • Overlapping: {mcl.get('overlapping', False)} → {leaf.get('overlapping', True)}")
        
        if pd.notna(mcl.get('mean_go_jaccard')) and pd.notna(leaf.get('mean_go_jaccard')):
            print(f"  • GO Jaccard: {mcl.get('mean_go_jaccard', 'N/A'):.6f} → {leaf.get('mean_go_jaccard', 'N/A'):.6f}")
            jaccard_improvement, _ = calculate_improvement(mcl.get('mean_go_jaccard'), leaf.get('mean_go_jaccard'), 'mean_go_jaccard', True)
            if jaccard_improvement is not None:
                print(f"    Improvement: {jaccard_improvement:+.2f}%")
    
    print()
    print("=" * 100)

def main():
    """Main comparison analysis."""
    print("Loading comparison data...")
    df_string, df_gavin = load_comparison_data()
    
    print("Extracting MCL vs LEAF-PPI comparisons...")
    string_data = extract_mcl_vs_leaf(df_string, 'STRING')
    gavin_data = extract_mcl_vs_leaf(df_gavin, 'Gavin')
    
    data = [string_data, gavin_data]
    
    print("Creating comparison tables...")
    comparison_df = create_comparison_table(data)
    summary_df = create_summary_statistics(data)
    
    # Save to CSV
    output_file = 'mcl_vs_leaf_comparison.csv'
    comparison_df.to_csv(output_file, index=False)
    print(f"\n✓ Comparison table saved to: {output_file}")
    
    summary_file = 'mcl_vs_leaf_summary.csv'
    summary_df.to_csv(summary_file, index=False)
    print(f"✓ Summary statistics saved to: {summary_file}")
    
    # Print report
    print_comparison_report(data, comparison_df, summary_df)
    
    # Additional analysis: show parameter optimization
    print("\nLEAF-PPI OPTIMIZATION PARAMETERS")
    print("-" * 100)
    for dataset_data in data:
        if dataset_data is None:
            continue
        
        dataset_name = dataset_data['dataset']
        leaf_params = json.loads(dataset_data['leaf'].get('parameters', '{}'))
        
        print(f"\n{dataset_name} Dataset:")
        print(f"  • Alpha (permanence weight): {leaf_params.get('alpha', 'N/A'):.4f}")
        print(f"  • Overlap threshold (tau): {leaf_params.get('overlap_tau', 'N/A'):.4f}")
        print(f"  • Transfer threshold: {leaf_params.get('transfer_tau', 'N/A'):.4f}")
        print(f"  • LEA evaluations: {leaf_params.get('lea_evaluations', 'N/A')}")
    
    print("\n" + "=" * 100)
    print("Comparison complete!")

if __name__ == '__main__':
    main()

