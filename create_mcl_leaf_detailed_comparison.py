#!/usr/bin/env python3
"""
Create detailed MCL vs LEAF-PPI comparison with additional metrics.

This script:
1. Loads actual cluster data from outputs
2. Computes detailed metrics (intra-density, inter-density, mean FD)
3. Creates comprehensive comparison tables
4. Generates LaTeX table for manuscript
"""

import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
import numpy as np
import networkx as nx
import json
from typing import Dict, Set

from compare_methods import load_string_dataset, load_gavin_dataset
from src.evaluation import (
    calculate_intra_density, calculate_inter_density, 
    calculate_conductance, calculate_overlapping_modularity,
    calculate_mean_fd_per_cluster, calculate_go_jaccard_similarity
)
from src.mcl_clustering import MCLClustering
from src.lea.optimize import optimize_communities
from src.go_tfidf import GOTFIDF
from src.permanence import calculate_permanence_all_proteins

def load_clusters_from_outputs(dataset_name, output_dir):
    """Load clusters from output CSV files."""
    mcl_file = f"{output_dir}/clusters_initial_mcl.csv"
    lea_file = f"{output_dir}/clusters_optimized_lea.csv"
    
    mcl_clusters = {}
    lea_clusters = {}
    
    # Load MCL clusters
    if os.path.exists(mcl_file):
        df_mcl = pd.read_csv(mcl_file)
        for _, row in df_mcl.iterrows():
            cid = row['cluster_id']
            pid = row['protein_id']
            if cid not in mcl_clusters:
                mcl_clusters[cid] = set()
            mcl_clusters[cid].add(pid)
    
    # Load LEAF-PPI clusters
    if os.path.exists(lea_file):
        df_lea = pd.read_csv(lea_file)
        for _, row in df_lea.iterrows():
            cid = row['cluster_id']
            pid = row['protein_id']
            if cid not in lea_clusters:
                lea_clusters[cid] = set()
            lea_clusters[cid].add(pid)
    
    return mcl_clusters, lea_clusters

def compute_detailed_metrics(clusters, graph, protein_go_terms=None, go_tfidf=None):
    """Compute detailed evaluation metrics."""
    metrics = {}
    
    # Structural metrics
    intra_densities = []
    inter_densities = []
    conductances = []
    
    for cluster_id, cluster in clusters.items():
        if len(cluster) == 0:
            continue
        
        intra_densities.append(calculate_intra_density(cluster, graph))
        inter_densities.append(calculate_inter_density(cluster, clusters, graph))
        conductances.append(calculate_conductance(cluster, graph))
    
    metrics['intra_density_mean'] = np.mean(intra_densities) if intra_densities else 0.0
    metrics['inter_density_mean'] = np.mean(inter_densities) if inter_densities else 0.0
    metrics['conductance_mean'] = np.mean(conductances) if conductances else 1.0
    
    # Modularity
    metrics['modularity'] = calculate_overlapping_modularity(clusters, graph)
    
    # Biological metrics
    if protein_go_terms and go_tfidf:
        metrics['mean_fd_per_cluster'] = calculate_mean_fd_per_cluster(
            clusters, protein_go_terms, go_tfidf
        )
        metrics['mean_go_jaccard'] = calculate_go_jaccard_similarity(
            clusters, protein_go_terms, None
        )
    else:
        metrics['mean_fd_per_cluster'] = None
        metrics['mean_go_jaccard'] = None
    
    # Cluster statistics
    cluster_sizes = [len(c) for c in clusters.values() if len(c) > 0]
    metrics['num_clusters'] = len(cluster_sizes)
    metrics['mean_cluster_size'] = np.mean(cluster_sizes) if cluster_sizes else 0.0
    metrics['max_cluster_size'] = max(cluster_sizes) if cluster_sizes else 0
    metrics['min_cluster_size'] = min(cluster_sizes) if cluster_sizes else 0
    
    # Overlapping statistics
    protein_cluster_count = {}
    for cluster in clusters.values():
        for protein in cluster:
            protein_cluster_count[protein] = protein_cluster_count.get(protein, 0) + 1
    
    overlapping_proteins = sum(1 for count in protein_cluster_count.values() if count > 1)
    metrics['overlapping_proteins'] = overlapping_proteins
    metrics['overlapping_percentage'] = (overlapping_proteins / len(protein_cluster_count) * 100) if protein_cluster_count else 0.0
    metrics['mean_clusters_per_protein'] = np.mean(list(protein_cluster_count.values())) if protein_cluster_count else 1.0
    
    return metrics

def create_detailed_comparison():
    """Create detailed comparison from actual cluster data."""
    results = []
    
    # Process STRING dataset
    print("Processing STRING dataset...")
    try:
        graph_str, lea_data_str = load_string_dataset()
        protein_go_terms_str = lea_data_str.get('protein_go_terms', {})
        go_tfidf_str = lea_data_str.get('go_tfidf')
        
        # Get MCL clusters (initial)
        mcl_clusters_str = lea_data_str.get('initial_clusters', {})
        
        # Get LEAF-PPI clusters (try to load from outputs, or use initial)
        lea_clusters_str = mcl_clusters_str  # Default to MCL if outputs not available
        if os.path.exists('outputs/clusters_optimized_lea.csv'):
            df_lea = pd.read_csv('outputs/clusters_optimized_lea.csv')
            lea_clusters_str = {}
            for _, row in df_lea.iterrows():
                cid = row['cluster_id']
                pid = row['protein_id']
                if cid not in lea_clusters_str:
                    lea_clusters_str[cid] = set()
                lea_clusters_str[cid].add(pid)
        
        mcl_metrics_str = compute_detailed_metrics(mcl_clusters_str, graph_str, protein_go_terms_str, go_tfidf_str)
        lea_metrics_str = compute_detailed_metrics(lea_clusters_str, graph_str, protein_go_terms_str, go_tfidf_str)
        
        results.append({
            'dataset': 'STRING',
            'method': 'MCL',
            **mcl_metrics_str
        })
        results.append({
            'dataset': 'STRING',
            'method': 'LEAF-PPI',
            **lea_metrics_str
        })
    except Exception as e:
        print(f"Error processing STRING: {e}")
    
    # Process Gavin dataset
    print("Processing Gavin dataset...")
    try:
        graph_gav, lea_data_gav = load_gavin_dataset()
        protein_go_terms_gav = lea_data_gav.get('protein_go_terms', {})
        go_tfidf_gav = lea_data_gav.get('go_tfidf')
        
        # Get MCL clusters
        mcl_clusters_gav = lea_data_gav.get('initial_clusters', {})
        
        # Get LEAF-PPI clusters
        lea_clusters_gav = mcl_clusters_gav
        if os.path.exists('outputs_gavin/clusters_optimized_lea.csv'):
            df_lea = pd.read_csv('outputs_gavin/clusters_optimized_lea.csv')
            lea_clusters_gav = {}
            for _, row in df_lea.iterrows():
                cid = row['cluster_id']
                pid = row['protein_id']
                if cid not in lea_clusters_gav:
                    lea_clusters_gav[cid] = set()
                lea_clusters_gav[cid].add(pid)
        
        mcl_metrics_gav = compute_detailed_metrics(mcl_clusters_gav, graph_gav, protein_go_terms_gav, go_tfidf_gav)
        lea_metrics_gav = compute_detailed_metrics(lea_clusters_gav, graph_gav, protein_go_terms_gav, go_tfidf_gav)
        
        results.append({
            'dataset': 'Gavin',
            'method': 'MCL',
            **mcl_metrics_gav
        })
        results.append({
            'dataset': 'Gavin',
            'method': 'LEAF-PPI',
            **lea_metrics_gav
        })
    except Exception as e:
        print(f"Error processing Gavin: {e}")
    
    df = pd.DataFrame(results)
    return df

def create_latex_table(df):
    """Create LaTeX table for manuscript."""
    latex_rows = []
    
    metrics_to_show = [
        ('modularity', 'Modularity', 4),
        ('intra_density_mean', 'Intra-Density', 4),
        ('inter_density_mean', 'Inter-Density', 4),
        ('conductance_mean', 'Conductance', 4),
        ('mean_fd_per_cluster', 'Mean FD', 4),
        ('mean_go_jaccard', 'GO Jaccard', 4),
        ('num_clusters', 'Communities', 0),
        ('overlapping_proteins', 'Overlapping Proteins', 0),
        ('overlapping_percentage', 'Overlap \%', 2),
    ]
    
    for dataset in ['STRING', 'Gavin']:
        df_dataset = df[df['dataset'] == dataset]
        mcl_row = df_dataset[df_dataset['method'] == 'MCL'].iloc[0]
        lea_row = df_dataset[df_dataset['method'] == 'LEAF-PPI'].iloc[0]
        
        for metric_key, metric_name, decimals in metrics_to_show:
            mcl_val = mcl_row.get(metric_key)
            lea_val = lea_row.get(metric_key)
            
            if pd.isna(mcl_val) or pd.isna(lea_val):
                continue
            
            if decimals == 0:
                mcl_str = f"{int(mcl_val)}"
                lea_str = f"{int(lea_val)}"
            else:
                mcl_str = f"{mcl_val:.{decimals}f}"
                lea_str = f"{lea_val:.{decimals}f}"
            
            # Calculate improvement
            if mcl_val != 0:
                improvement = ((lea_val - mcl_val) / abs(mcl_val)) * 100
                improvement_str = f"{improvement:+.1f}\\%"
            else:
                improvement_str = "---"
            
            latex_rows.append({
                'Dataset': dataset,
                'Metric': metric_name,
                'MCL': mcl_str,
                'LEAF-PPI': lea_str,
                'Improvement': improvement_str
            })
    
    return pd.DataFrame(latex_rows)

def main():
    """Main function."""
    print("=" * 80)
    print("Creating Detailed MCL vs LEAF-PPI Comparison")
    print("=" * 80)
    
    df = create_detailed_comparison()
    
    # Save detailed comparison
    output_file = 'mcl_vs_leaf_detailed.csv'
    df.to_csv(output_file, index=False)
    print(f"\n✓ Detailed comparison saved to: {output_file}")
    
    # Create LaTeX table
    latex_df = create_latex_table(df)
    latex_file = 'mcl_vs_leaf_latex_table.csv'
    latex_df.to_csv(latex_file, index=False)
    print(f"✓ LaTeX table data saved to: {latex_file}")
    
    # Print summary
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY")
    print("=" * 80)
    print(df.to_string())
    
    print("\n" + "=" * 80)
    print("KEY IMPROVEMENTS")
    print("=" * 80)
    
    for dataset in ['STRING', 'Gavin']:
        df_dataset = df[df['dataset'] == dataset]
        mcl_row = df_dataset[df_dataset['method'] == 'MCL'].iloc[0]
        lea_row = df_dataset[df_dataset['method'] == 'LEAF-PPI'].iloc[0]
        
        print(f"\n{dataset} Dataset:")
        print(f"  Modularity: {mcl_row['modularity']:.6f} → {lea_row['modularity']:.6f}")
        if mcl_row['modularity'] != 0:
            mod_imp = ((lea_row['modularity'] - mcl_row['modularity']) / abs(mcl_row['modularity'])) * 100
            print(f"    Improvement: {mod_imp:+.2f}%")
        
        print(f"  Overlapping proteins: {mcl_row['overlapping_proteins']} → {lea_row['overlapping_proteins']}")
        print(f"  Overlap percentage: {mcl_row['overlapping_percentage']:.2f}% → {lea_row['overlapping_percentage']:.2f}%")
        
        if pd.notna(mcl_row.get('mean_fd_per_cluster')) and pd.notna(lea_row.get('mean_fd_per_cluster')):
            print(f"  Mean FD: {mcl_row['mean_fd_per_cluster']:.6f} → {lea_row['mean_fd_per_cluster']:.6f}")

if __name__ == '__main__':
    main()

