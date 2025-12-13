#!/usr/bin/env python3
"""
Generate updated CSV result files with normalized metrics and external GO evaluation.

This script re-runs all experiments with:
- Normalized permanence [-1, 1]
- Normalized FD [-1, 1]
- MCL filtering (min_cluster_size=10)
- External GO Jaccard evaluation
- Fixed random seeds for reproducibility
"""

import logging
import sys
import os
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent))

from compare_methods import load_string_dataset, load_gavin_dataset
from src.community_comparison import compare_all_methods

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Generate updated CSV files for both datasets."""
    
    logger.info("=" * 80)
    logger.info("Generating Updated Results with Normalized Metrics")
    logger.info("=" * 80)
    logger.info("\nKey Updates:")
    logger.info("  - Permanence normalized to [-1, 1]")
    logger.info("  - Functional Dependency normalized to [-1, 1]")
    logger.info("  - MCL filtering: min_cluster_size=10")
    logger.info("  - External GO evaluation: mean_go_jaccard")
    logger.info("  - Fixed random seed: 42")
    logger.info("=" * 80)
    
    random_seed = 42
    lea_evaluations = 500
    
    all_results = []
    
    # Dataset 1: STRING
    logger.info("\n" + "=" * 80)
    logger.info("DATASET 1: STRING")
    logger.info("=" * 80)
    try:
        graph_str, lea_data_str = load_string_dataset()
        protein_go_terms_str = lea_data_str.get('protein_go_terms', {}) if lea_data_str else {}
        
        results_str = compare_all_methods(
            graph_str,
            'STRING',
            ground_truth=None,
            lea_data=lea_data_str,
            lea_evaluations=lea_evaluations,
            random_seed=random_seed,
            protein_go_terms=protein_go_terms_str
        )
        all_results.append(results_str)
        logger.info(f"✓ STRING: {len(results_str)} methods completed")
    except Exception as e:
        logger.error(f"✗ STRING dataset failed: {e}", exc_info=True)
    
    # Dataset 2: Gavin
    logger.info("\n" + "=" * 80)
    logger.info("DATASET 2: GAVIN")
    logger.info("=" * 80)
    try:
        graph_gav, lea_data_gav = load_gavin_dataset()
        protein_go_terms_gav = lea_data_gav.get('protein_go_terms', {}) if lea_data_gav else {}
        
        results_gav = compare_all_methods(
            graph_gav,
            'Gavin',
            ground_truth=None,
            lea_data=lea_data_gav,
            lea_evaluations=lea_evaluations,
            random_seed=random_seed,
            protein_go_terms=protein_go_terms_gav
        )
        all_results.append(results_gav)
        logger.info(f"✓ Gavin: {len(results_gav)} methods completed")
    except Exception as e:
        logger.error(f"✗ Gavin dataset failed: {e}", exc_info=True)
    
    # Save results
    if all_results:
        import pandas as pd
        
        # Combine all results
        df = pd.concat(all_results, ignore_index=True)
        
        # Save combined results
        output_file = 'community_detection_comparison.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"\n{'=' * 80}")
        logger.info(f"✓ Combined results saved to: {output_file}")
        
        # Save separate files per dataset
        for dataset_name in ['STRING', 'Gavin']:
            df_dataset = df[df['dataset'] == dataset_name]
            if len(df_dataset) > 0:
                output_file_dataset = f'results_{dataset_name.lower()}_updated.csv'
                df_dataset.to_csv(output_file_dataset, index=False)
                logger.info(f"✓ {dataset_name} results saved to: {output_file_dataset}")
        
        logger.info(f"\nTotal comparisons: {len(df)}")
        logger.info(f"{'=' * 80}")
        
        # Print summary
        print("\n" + "=" * 80)
        print("SUMMARY")
        print("=" * 80)
        print(df.to_string())
        print("=" * 80)
        
    else:
        logger.error("✗ No results generated!")
        sys.exit(1)


if __name__ == '__main__':
    main()

