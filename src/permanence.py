"""
Permanence calculation for community detection.
Implements Eq.1 from main.docx.

Permanence measures how well a protein belongs to its community based on
internal vs external connections.
"""

import logging
from typing import Dict, Set
import networkx as nx

logger = logging.getLogger(__name__)


def calculate_permanence(protein: str, cluster: Set[str], graph: nx.Graph,
                         all_clusters: Dict[int, Set[str]] = None) -> float:
    """
    Calculate permanence for a protein in a cluster (Eq.1).
    
    Permanence(p, c) = (I(p) / E_max(p)) - (1 - C_in(p))
    
    Where:
    - I(p): internal connections (edges to proteins in same cluster)
    - E_max(p): maximum external connections to any single cluster
    - C_in(p): clustering coefficient of internal neighbors
    
    IMPORTANT: Permanence is normalized to [-1, 1] range to ensure
    theoretical consistency with literature definitions.
    
    Args:
        protein: Protein ID
        cluster: Set of protein IDs in the cluster
        graph: NetworkX graph
        all_clusters: Dict mapping cluster_id to set of proteins (for E_max calculation)
        
    Returns:
        Permanence score in range [-1, 1]
    """
    if protein not in graph:
        return 0.0
    
    neighbors = set(graph.neighbors(protein))
    if not neighbors:
        return 0.0
    
    # Internal connections: neighbors in the same cluster
    internal_neighbors = neighbors & cluster
    I_p = len(internal_neighbors)
    
    # External connections: neighbors not in this cluster
    external_neighbors = neighbors - cluster
    
    if not external_neighbors:
        # No external connections - fully internal
        if I_p > 0:
            # Normalize to ensure range [-1, 1]
            # With no external connections, permanence should be positive
            return 1.0
        return 0.0
    
    # E_max: maximum external connections to any single OTHER cluster
    # This requires knowledge of all clusters
    E_max_p = 0
    if all_clusters is not None:
        for cluster_id, other_cluster in all_clusters.items():
            if other_cluster == cluster:
                continue  # Skip the current cluster
            connections = len(neighbors & other_cluster)
            if connections > E_max_p:
                E_max_p = connections
    else:
        # Fallback: use total external neighbors (upper bound)
        E_max_p = len(external_neighbors)
    
    # Clustering coefficient of internal neighbors
    # C_in(p) = fraction of edges between internal neighbors
    if I_p < 2:
        C_in_p = 0.0
    else:
        internal_edges = 0
        internal_list = list(internal_neighbors)
        for i, n1 in enumerate(internal_list):
            for n2 in internal_list[i+1:]:
                if graph.has_edge(n1, n2):
                    internal_edges += 1
        max_possible = I_p * (I_p - 1) / 2
        C_in_p = internal_edges / max_possible if max_possible > 0 else 0.0
    
    # Calculate permanence (Eq.1)
    if E_max_p == 0:
        # No external connections to any cluster
        permanence = (I_p / max(1, len(neighbors))) - (1 - C_in_p)
    else:
        permanence = (I_p / E_max_p) - (1 - C_in_p)
    
    # Normalize permanence to [-1, 1] range
    # Theoretical bounds: 
    # - Lower bound: when I_p = 0 and C_in_p = 0, permanence = -1
    # - Upper bound: when I_p >> E_max_p and C_in_p = 1, permanence approaches 1
    # Use linear clipping to ensure strict bounds (more conservative and matches theory)
    permanence_normalized = max(-1.0, min(1.0, permanence))
    
    return permanence_normalized


def calculate_permanence_all_proteins(clusters: Dict[int, Set[str]], 
                                      graph: nx.Graph) -> Dict[str, Dict[int, float]]:
    """
    Calculate permanence for all proteins in all clusters.
    
    Args:
        clusters: Dict mapping cluster_id to set of protein IDs
        graph: NetworkX graph
        
    Returns:
        Dict mapping protein_id to dict of cluster_id -> permanence score
    """
    permanence_scores = {}
    
    # Build reverse mapping: protein -> clusters it belongs to
    protein_clusters = {}
    for cluster_id, proteins in clusters.items():
        for protein in proteins:
            if protein not in protein_clusters:
                protein_clusters[protein] = []
            protein_clusters[protein].append(cluster_id)
    
    # Calculate permanence for each protein in each cluster
    # Pass all_clusters to properly compute E_max
    for protein, cluster_ids in protein_clusters.items():
        permanence_scores[protein] = {}
        for cluster_id in cluster_ids:
            cluster = clusters[cluster_id]
            perm = calculate_permanence(protein, cluster, graph, all_clusters=clusters)
            permanence_scores[protein][cluster_id] = perm
    
    return permanence_scores

