"""
Fitness function for LEA optimization of community membership.
Maximizes community quality metrics while minimizing fragmentation.
"""

import logging
from typing import Dict, Set, Callable
import networkx as nx
import numpy as np

logger = logging.getLogger(__name__)


class MembershipFitness:
    """
    Fitness function for optimizing community membership parameters.
    """
    
    def __init__(self, graph: nx.Graph,
                 initial_clusters: Dict[int, Set[str]],
                 protein_go_terms: Dict[str, Set[str]],
                 go_tfidf,
                 permanence_scores: Dict[str, Dict[int, float]]):
        """
        Initialize fitness function.
        
        Args:
            graph: NetworkX PPI graph
            initial_clusters: Initial MCL clusters
            protein_go_terms: Dict mapping protein ID to GO terms
            go_tfidf: GOTFIDF instance
            permanence_scores: Pre-computed permanence scores
        """
        self.graph = graph
        self.initial_clusters = initial_clusters
        self.protein_go_terms = protein_go_terms
        self.go_tfidf = go_tfidf
        self.permanence_scores = permanence_scores
        
        # Cache for performance
        self._membership_cache = {}
    
    def compute_fitness(self, solution: np.ndarray, 
                       lambda_inter: float = 1.0,
                       lambda_fragment: float = 0.5) -> float:
        """
        Compute fitness for a solution vector.
        
        Solution vector components (example):
        [0] alpha: Membership weight (0-1)
        [1] overlap_tau: Overlap threshold (0-1)
        [2] transfer_tau: Transfer threshold (0-1)
        
        Fitness increases with:
        - Average membership scores
        - Intra-cluster cohesion
        - GO coherence
        
        Fitness decreases with:
        - Inter-cluster coupling
        - Fragmentation (too many small clusters)
        
        Args:
            solution: Decision variables vector
            lambda_inter: Weight for inter-cluster penalty
            lambda_fragment: Weight for fragmentation penalty
            
        Returns:
            Fitness value (to maximize, so we'll negate for minimization)
        """
        # Extract parameters from solution
        alpha = np.clip(solution[0], 0.0, 1.0)
        overlap_tau = np.clip(solution[1], 0.0, 1.0)
        transfer_tau = np.clip(solution[2], 0.0, 1.0)
        
        # Apply overlap reassignment with these parameters
        from src.membership_overlap import apply_overlap_reassignment
        
        try:
            optimized_clusters = apply_overlap_reassignment(
                self.initial_clusters,
                self.graph,
                self.protein_go_terms,
                self.go_tfidf,
                self.permanence_scores,
                alpha=alpha,
                overlap_tau=overlap_tau,
                transfer_tau=transfer_tau
            )
        except Exception as e:
            logger.warning(f"Error in overlap reassignment: {e}")
            return -1e6  # Very bad fitness
        
        # Compute fitness components
        membership_sum = 0.0
        membership_count = 0
        
        intra_cohesion = 0.0
        inter_coupling = 0.0
        
        go_coherence_sum = 0.0
        go_coherence_count = 0
        
        # Calculate metrics for each cluster
        for cluster_id, cluster in optimized_clusters.items():
            if len(cluster) == 0:
                continue
            
            # Average membership in cluster
            cluster_membership_sum = 0.0
            cluster_size = len(cluster)
            
            # Intra-cluster edges
            intra_edges = 0
            cluster_list = list(cluster)
            for i, p1 in enumerate(cluster_list):
                for p2 in cluster_list[i+1:]:
                    if self.graph.has_edge(p1, p2):
                        intra_edges += 1
                
                # Calculate membership for this protein
                from src.membership_overlap import calculate_membership
                memb = calculate_membership(
                    p1, cluster, cluster_id, self.graph,
                    self.protein_go_terms, self.go_tfidf,
                    self.permanence_scores, alpha
                )
                cluster_membership_sum += memb
                membership_count += 1
            
            avg_membership = cluster_membership_sum / cluster_size if cluster_size > 0 else 0.0
            membership_sum += avg_membership * cluster_size
            
            # Intra-cohesion: fraction of possible edges that exist
            max_possible_edges = cluster_size * (cluster_size - 1) / 2
            if max_possible_edges > 0:
                cohesion = intra_edges / max_possible_edges
                intra_cohesion += cohesion * cluster_size
            
            # Inter-cluster coupling: edges from cluster to other clusters
            inter_edges = 0
            for protein in cluster:
                neighbors = set(self.graph.neighbors(protein))
                inter_edges += len(neighbors - cluster)
            
            inter_coupling += inter_edges
            
            # GO coherence: average functional dependency in cluster
            from src.membership_overlap import calculate_functional_dependency
            cluster_fd_sum = 0.0
            for protein in cluster:
                fd = calculate_functional_dependency(
                    protein, cluster, self.protein_go_terms,
                    self.go_tfidf, cluster_id
                )
                cluster_fd_sum += fd
            
            avg_fd = cluster_fd_sum / cluster_size if cluster_size > 0 else 0.0
            go_coherence_sum += avg_fd * cluster_size
            go_coherence_count += cluster_size
        
        # Normalize metrics
        total_proteins = sum(len(c) for c in optimized_clusters.values())
        if total_proteins == 0:
            return -1e6
        
        avg_membership = membership_sum / total_proteins if total_proteins > 0 else 0.0
        avg_intra_cohesion = intra_cohesion / total_proteins if total_proteins > 0 else 0.0
        avg_go_coherence = go_coherence_sum / total_proteins if total_proteins > 0 else 0.0
        
        # Normalize inter-coupling (per edge)
        total_edges = self.graph.number_of_edges()
        normalized_inter_coupling = inter_coupling / total_edges if total_edges > 0 else 0.0
        
        # Fragmentation penalty: penalize too many small clusters
        num_clusters = len([c for c in optimized_clusters.values() if len(c) > 0])
        num_singletons = len([c for c in optimized_clusters.values() if len(c) == 1])
        fragmentation = num_singletons / num_clusters if num_clusters > 0 else 1.0
        
        # Compute fitness (to maximize)
        fitness = (avg_membership + 
                  avg_intra_cohesion + 
                  avg_go_coherence -
                  lambda_inter * normalized_inter_coupling -
                  lambda_fragment * fragmentation)
        
        return fitness
    
    def create_fitness_function(self, lambda_inter: float = 1.0,
                                lambda_fragment: float = 0.5) -> Callable:
        """
        Create a fitness function compatible with LEA (minimization).
        
        Args:
            lambda_inter: Weight for inter-cluster penalty
            lambda_fragment: Weight for fragmentation penalty
            
        Returns:
            Fitness function that takes solution vector and returns scalar (to minimize)
        """
        def fitness_func(solution: np.ndarray) -> float:
            fitness = self.compute_fitness(solution, lambda_inter, lambda_fragment)
            # Negate for minimization (LEA minimizes)
            return -fitness
        
        return fitness_func

