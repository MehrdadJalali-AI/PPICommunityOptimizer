#!/usr/bin/env python3
"""
Main script to compare community detection methods across datasets.
Generates community_detection_comparison.csv
"""

import sys
import os
import logging
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.community_comparison import compare_all_methods
from src.string_loader import STRINGLoader
from src.gavin_loader import GavinLoader
from src.go_loader import GOLoader
from src.mcl_clustering import MCLClustering
from src.go_tfidf import GOTFIDF
from src.permanence import calculate_permanence_all_proteins
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_string_dataset(taxid: int = 4932, threshold: int = 700):
    """Load STRING dataset."""
    logger.info("Loading STRING dataset...")
    loader = STRINGLoader(taxid, cache_dir='cache', threshold=threshold)
    graph, aliases = loader.load_from_download()
    
    # Get initial clusters for LEA
    mcl = MCLClustering(inflation=2.0)
    initial_clusters = mcl.cluster(graph)
    
    # Try to load GO annotations
    go_loader = GOLoader('cache')
    gaf_file = f'cache/goa_{taxid}.gaf.gz'
    protein_go_terms = {}
    go_tfidf = None
    permanence_scores = {}
    
    if os.path.exists(gaf_file):
        protein_go_terms = go_loader.load_from_gaf(gaf_file, taxid)
        go_tfidf = GOTFIDF(initial_clusters, protein_go_terms)
        permanence_scores = calculate_permanence_all_proteins(initial_clusters, graph)
    
    lea_data = {
        'initial_clusters': initial_clusters,
        'protein_go_terms': protein_go_terms,
        'go_tfidf': go_tfidf,
        'permanence_scores': permanence_scores
    }
    
    return graph, lea_data


def load_gavin_dataset(ppi_file: str = 'gavin2006_socioaffinities_rescaled.txt',
                      go_file: str = 'GO.txt'):
    """Load Gavin dataset."""
    logger.info("Loading Gavin dataset...")
    loader = GavinLoader(normalize=True)
    graph = loader.load(ppi_file)
    
    # Get initial clusters for LEA
    mcl = MCLClustering(inflation=2.0)
    initial_clusters = mcl.cluster(graph)
    
    # Load GO annotations
    go_loader = GOLoader('cache')
    protein_go_terms = {}
    go_tfidf = None
    permanence_scores = {}
    
    if os.path.exists(go_file):
        protein_go_terms = go_loader.load_from_gaf(go_file, taxid=559292, use_symbol=True)
        go_tfidf = GOTFIDF(initial_clusters, protein_go_terms)
        permanence_scores = calculate_permanence_all_proteins(initial_clusters, graph)
    
    lea_data = {
        'initial_clusters': initial_clusters,
        'protein_go_terms': protein_go_terms,
        'go_tfidf': go_tfidf,
        'permanence_scores': permanence_scores
    }
    
    return graph, lea_data


def main():
    """Run comparison across all methods and datasets."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Compare community detection methods')
    parser.add_argument('--lea-evaluations', type=int, default=500,
                       help='Number of LEA evaluations (default: 500, use lower for faster testing)')
    parser.add_argument('--skip-lea', action='store_true',
                       help='Skip LEA method (for faster testing)')
    args = parser.parse_args()
    
    logger.info("=" * 60)
    logger.info("Community Detection Method Comparison")
    logger.info("=" * 60)
    
    all_results = []
    random_seed = 42
    
    # Dataset 1: STRING
    try:
        logger.info("\n" + "=" * 60)
        logger.info("DATASET 1: STRING")
        logger.info("=" * 60)
        graph_str, lea_data_str = load_string_dataset()
        if args.skip_lea:
            lea_data_str = None
        results_str = compare_all_methods(
            graph_str,
            'STRING',
            ground_truth=None,  # No ground truth available
            lea_data=lea_data_str,
            lea_evaluations=args.lea_evaluations,
            random_seed=random_seed
        )
        all_results.append(results_str)
        logger.info(f"STRING: {len(results_str)} methods completed")
    except Exception as e:
        logger.error(f"STRING dataset failed: {e}", exc_info=True)
    
    # Dataset 2: Gavin
    try:
        logger.info("\n" + "=" * 60)
        logger.info("DATASET 2: GAVIN")
        logger.info("=" * 60)
        graph_gav, lea_data_gav = load_gavin_dataset()
        if args.skip_lea:
            lea_data_gav = None
        results_gav = compare_all_methods(
            graph_gav,
            'Gavin',
            ground_truth=None,  # No ground truth available
            lea_data=lea_data_gav,
            lea_evaluations=args.lea_evaluations,
            random_seed=random_seed
        )
        all_results.append(results_gav)
        logger.info(f"Gavin: {len(results_gav)} methods completed")
    except Exception as e:
        logger.error(f"Gavin dataset failed: {e}", exc_info=True)
    
    # Combine results
    if all_results:
        df = pd.concat(all_results, ignore_index=True)
        
        # Save to CSV
        output_file = 'community_detection_comparison.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Results saved to: {output_file}")
        logger.info(f"Total comparisons: {len(df)}")
        logger.info(f"{'=' * 60}")
        
        # Print summary
        print("\n" + "=" * 60)
        print("COMPARISON SUMMARY")
        print("=" * 60)
        print(df.to_string())
        print("=" * 60)
        
        # Print best methods per metric
        print("\nBest Methods by Metric:")
        print("-" * 60)
        for metric in ['modularity', 'conductance']:
            if metric in df.columns:
                best = df.loc[df[metric].idxmax() if df[metric].max() > 0 else df[metric].idxmin()]
                print(f"{metric}: {best['method']} ({best['dataset']}) = {best[metric]:.4f}")
        
        if 'runtime_sec' in df.columns:
            fastest = df.loc[df['runtime_sec'].idxmin()]
            print(f"Fastest: {fastest['method']} ({fastest['dataset']}) = {fastest['runtime_sec']:.2f}s")
        
    else:
        logger.error("No results generated!")


if __name__ == '__main__':
    main()

