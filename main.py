#!/usr/bin/env python3
"""
Main CLI for community detection pipeline.

Usage:
    python main.py --taxid 4932 --string-mode download --threshold 700 --go-source goa --outdir outputs/
"""

import argparse
import logging
import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.string_loader import STRINGLoader
from src.gavin_loader import GavinLoader
from src.go_loader import GOLoader
from src.mcl_clustering import MCLClustering
from src.go_tfidf import GOTFIDF
from src.permanence import calculate_permanence_all_proteins
from src.membership_overlap import apply_overlap_reassignment
from src.lea.optimize import optimize_communities
from src.evaluation import evaluate_clusters
from src.outputs import (
    save_initial_clusters, save_go_term_importance,
    save_protein_membership, save_optimized_clusters,
    save_overlap_summary, save_evaluation_results
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Community Detection Pipeline with LEA Optimization'
    )
    
    # Data source selection
    parser.add_argument('--mode', choices=['string', 'gavin'], default='string',
                       help='Data source mode: string (STRING DB) or gavin (Gavin PPI)')
    
    # STRING-specific arguments
    parser.add_argument('--taxid', type=int, default=None,
                       help='NCBI taxonomy ID (e.g., 4932 for S. cerevisiae, required for STRING mode)')
    parser.add_argument('--string-mode', choices=['download', 'api'], default='download',
                       help='STRING data loading mode')
    parser.add_argument('--threshold', type=int, default=700,
                       help='STRING combined score threshold (0-1000)')
    
    # Gavin-specific arguments
    parser.add_argument('--ppi', type=str, default=None,
                       help='Path to PPI file (Gavin format: protein1\\tprotein2\\tweight)')
    
    # GO arguments
    parser.add_argument('--go-file', type=str, required=True,
                       help='Path to GO annotation file (GAF format)')
    parser.add_argument('--go-use-symbol', action='store_true',
                       help='Use DB_Object_Symbol instead of DB_Object_ID (for SGD GAF files)')
    parser.add_argument('--go-taxid', type=int, default=None,
                       help='Taxonomy ID to filter GO annotations (e.g., 559292 for yeast)')
    
    # Common arguments
    parser.add_argument('--outdir', type=str, default='outputs',
                       help='Output directory')
    parser.add_argument('--cache-dir', type=str, default='cache',
                       help='Cache directory for downloaded files')
    parser.add_argument('--mcl-inflation', type=float, default=2.0,
                       help='MCL inflation parameter')
    parser.add_argument('--alpha', type=float, default=0.5,
                       help='Membership weight parameter (0-1)')
    parser.add_argument('--overlap-tau', type=float, default=0.1,
                       help='Overlap threshold for adding proteins to additional clusters')
    parser.add_argument('--transfer-tau', type=float, default=0.0,
                       help='Transfer threshold')
    parser.add_argument('--lea-population', type=int, default=30,
                       help='LEA population size')
    parser.add_argument('--lea-evaluations', type=int, default=500,
                       help='LEA maximum function evaluations')
    parser.add_argument('--lambda-inter', type=float, default=1.0,
                       help='Weight for inter-cluster penalty')
    parser.add_argument('--lambda-fragment', type=float, default=0.5,
                       help='Weight for fragmentation penalty')
    parser.add_argument('--random-seed', type=int, default=None,
                       help='Random seed for reproducibility')
    parser.add_argument('--skip-lea', action='store_true',
                       help='Skip LEA optimization (use initial MCL clusters)')
    parser.add_argument('--gold-standard', type=str, default=None,
                       help='Path to gold standard clusters file (optional)')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    logger.info("=" * 60)
    logger.info("Community Detection Pipeline")
    logger.info("=" * 60)
    logger.info(f"Mode: {args.mode.upper()}")
    logger.info(f"Output directory: {args.outdir}")
    logger.info("=" * 60)
    
    # Step 1: Load PPI network
    logger.info("\n[Step 1] Loading PPI network...")
    
    if args.mode == 'gavin':
        if not args.ppi:
            logger.error("--ppi argument required for Gavin mode")
            sys.exit(1)
        
        gavin_loader = GavinLoader(normalize=True)
        graph = gavin_loader.load(args.ppi)
        protein_aliases = {}  # No aliases for Gavin format
        
        logger.info(f"Loaded Gavin network: {graph.number_of_nodes()} nodes, "
                   f"{graph.number_of_edges()} edges")
    else:  # STRING mode
        if not args.taxid:
            logger.error("--taxid argument required for STRING mode")
            sys.exit(1)
        
        string_loader = STRINGLoader(args.taxid, args.cache_dir, args.threshold)
        
        if args.string_mode == 'download':
            graph, protein_aliases = string_loader.load_from_download()
        else:
            graph, protein_aliases = string_loader.load_from_api()
        
        logger.info(f"Loaded STRING network: {graph.number_of_nodes()} nodes, "
                   f"{graph.number_of_edges()} edges")
    
    # Step 2: Load GO annotations
    logger.info("\n[Step 2] Loading GO annotations...")
    go_loader = GOLoader(args.cache_dir)
    
    if not os.path.exists(args.go_file):
        logger.error(f"GO file not found: {args.go_file}")
        sys.exit(1)
    
    protein_go_terms = go_loader.load_from_gaf(
        args.go_file, 
        taxid=args.go_taxid,
        use_symbol=args.go_use_symbol
    )
    
    logger.info(f"Loaded GO annotations for {len(protein_go_terms)} proteins")
    
    # Step 3: MCL clustering
    logger.info("\n[Step 3] Running MCL clustering...")
    mcl = MCLClustering(inflation=args.mcl_inflation)
    initial_clusters = mcl.cluster(graph)
    logger.info(f"MCL found {len(initial_clusters)} clusters")
    
    # Save initial clusters
    save_initial_clusters(
        initial_clusters,
        os.path.join(args.outdir, 'clusters_initial_mcl.csv')
    )
    
    # Step 4: Compute GO TF-IDF
    logger.info("\n[Step 4] Computing GO TF-IDF importance...")
    go_tfidf = GOTFIDF(initial_clusters, protein_go_terms)
    
    # Save GO term importance
    save_go_term_importance(
        go_tfidf,
        os.path.join(args.outdir, 'go_term_importance.csv')
    )
    
    # Step 5: Calculate permanence
    logger.info("\n[Step 5] Calculating permanence...")
    permanence_scores = calculate_permanence_all_proteins(initial_clusters, graph)
    
    # Step 6: Calculate membership and apply overlap (initial)
    logger.info("\n[Step 6] Applying initial overlap reassignment...")
    
    clusters_with_overlap = apply_overlap_reassignment(
        initial_clusters,
        graph,
        protein_go_terms,
        go_tfidf,
        permanence_scores,
        alpha=args.alpha,
        overlap_tau=args.overlap_tau,
        transfer_tau=args.transfer_tau
    )
    
    # Save protein membership
    save_protein_membership(
        clusters_with_overlap,
        graph,
        protein_go_terms,
        go_tfidf,
        permanence_scores,
        args.alpha,
        os.path.join(args.outdir, 'protein_membership.csv')
    )
    
    # Step 7: LEA optimization
    if not args.skip_lea:
        logger.info("\n[Step 7] Running LEA optimization...")
        best_solution, best_fitness, optimized_clusters = optimize_communities(
            graph,
            initial_clusters,
            protein_go_terms,
            go_tfidf,
            permanence_scores,
            population_size=args.lea_population,
            max_evaluations=args.lea_evaluations,
            lambda_inter=args.lambda_inter,
            lambda_fragment=args.lambda_fragment,
            random_seed=args.random_seed
        )
        
        optimized_alpha = best_solution[0]
    else:
        logger.info("\n[Step 7] Skipping LEA optimization (using initial clusters)")
        optimized_clusters = clusters_with_overlap
        optimized_alpha = args.alpha
    
    # Save optimized clusters
    save_optimized_clusters(
        optimized_clusters,
        permanence_scores,
        protein_go_terms,
        go_tfidf,
        graph,
        optimized_alpha,
        os.path.join(args.outdir, 'clusters_optimized_lea.csv')
    )
    
    # Save overlap summary
    save_overlap_summary(
        optimized_clusters,
        os.path.join(args.outdir, 'overlap_summary.csv')
    )
    
    # Step 8: Evaluation
    logger.info("\n[Step 8] Evaluating clusters...")
    
    gold_standard = None
    if args.gold_standard and os.path.exists(args.gold_standard):
        logger.info(f"Loading gold standard from {args.gold_standard}")
        import pandas as pd
        gold_df = pd.read_csv(args.gold_standard)
        gold_standard = {}
        for _, row in gold_df.iterrows():
            cid = row['cluster_id']
            pid = row['protein_id']
            if cid not in gold_standard:
                gold_standard[cid] = set()
            gold_standard[cid].add(pid)
    
    evaluation_df = evaluate_clusters(
        optimized_clusters,
        graph,
        protein_go_terms,
        go_tfidf,
        gold_standard=gold_standard
    )
    
    save_evaluation_results(
        evaluation_df,
        os.path.join(args.outdir, 'evaluation_results.csv')
    )
    
    logger.info("\n" + "=" * 60)
    logger.info("Pipeline completed successfully!")
    logger.info("=" * 60)
    logger.info(f"\nResults saved to: {args.outdir}/")
    logger.info("\nEvaluation Summary:")
    print(evaluation_df.to_string())


if __name__ == '__main__':
    main()

