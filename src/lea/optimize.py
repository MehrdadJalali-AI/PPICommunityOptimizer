"""
LEA optimization wrapper for community detection.
"""

import logging
from typing import Dict, Set, Tuple
import numpy as np
import networkx as nx

from .lotus_effect_algorithm import LotusEffectAlgorithm
from .fitness_membership import MembershipFitness

logger = logging.getLogger(__name__)


def optimize_communities(graph: nx.Graph,
                         initial_clusters: Dict[int, Set[str]],
                         protein_go_terms: Dict[str, Set[str]],
                         go_tfidf,
                         permanence_scores: Dict[str, Dict[int, float]],
                         population_size: int = 30,
                         max_evaluations: int = 1000,
                         lambda_inter: float = 1.0,
                         lambda_fragment: float = 0.5,
                         random_seed: int = None) -> Tuple[np.ndarray, float, Dict[int, Set[str]]]:
    """
    Optimize community membership using LEA.
    
    Decision variables:
    [0] alpha: Membership weight (0-1)
    [1] overlap_tau: Overlap threshold (0-1)
    [2] transfer_tau: Transfer threshold (0-1)
    
    Args:
        graph: NetworkX PPI graph
        initial_clusters: Initial MCL clusters
        protein_go_terms: Dict mapping protein ID to GO terms
        go_tfidf: GOTFIDF instance
        permanence_scores: Pre-computed permanence scores
        population_size: LEA population size
        max_evaluations: Maximum function evaluations
        lambda_inter: Weight for inter-cluster penalty
        lambda_fragment: Weight for fragmentation penalty
        random_seed: Random seed for reproducibility
        
    Returns:
        (best_solution, best_fitness, optimized_clusters)
    """
    logger.info("Starting LEA optimization for community membership...")
    
    # Create fitness function
    fitness_obj = MembershipFitness(
        graph, initial_clusters, protein_go_terms,
        go_tfidf, permanence_scores
    )
    fitness_func = fitness_obj.create_fitness_function(lambda_inter, lambda_fragment)
    
    # Define decision variables (3 dimensions)
    dimensions = 3
    lower_bounds = np.array([0.0, 0.0, 0.0])  # alpha, overlap_tau, transfer_tau
    upper_bounds = np.array([1.0, 1.0, 1.0])
    
    # Initialize and run LEA
    lea = LotusEffectAlgorithm(
        population_size=population_size,
        dimensions=dimensions,
        lower_bound=lower_bounds,
        upper_bound=upper_bounds,
        max_function_evaluations=max_evaluations,
        fitness_function=fitness_func,
        random_seed=random_seed
    )
    
    best_solution, best_fitness, evaluations = lea.optimize()
    
    logger.info(f"LEA optimization completed: {evaluations} evaluations")
    logger.info(f"Best solution: alpha={best_solution[0]:.3f}, "
               f"overlap_tau={best_solution[1]:.3f}, "
               f"transfer_tau={best_solution[2]:.3f}")
    logger.info(f"Best fitness: {best_fitness:.6f}")
    
    # Apply best solution to get optimized clusters
    alpha = np.clip(best_solution[0], 0.0, 1.0)
    overlap_tau = np.clip(best_solution[1], 0.0, 1.0)
    transfer_tau = np.clip(best_solution[2], 0.0, 1.0)
    
    from src.membership_overlap import apply_overlap_reassignment
    
    optimized_clusters = apply_overlap_reassignment(
        initial_clusters,
        graph,
        protein_go_terms,
        go_tfidf,
        permanence_scores,
        alpha=alpha,
        overlap_tau=overlap_tau,
        transfer_tau=transfer_tau
    )
    
    return best_solution, best_fitness, optimized_clusters

