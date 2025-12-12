"""
Membership and functional dependency calculations.
Implements Eq.2 (functional dependency fd) and Eq.4 (Membership).
Handles overlapping community assignment.
"""

import logging
from typing import Dict, Set, List, Tuple
import networkx as nx

logger = logging.getLogger(__name__)


def calculate_functional_dependency(protein: str, cluster: Set[str],
                                   protein_go_terms: Dict[str, Set[str]],
                                   go_tfidf: 'GOTFIDF',
                                   cluster_id: int) -> float:
    """
    Calculate functional dependency fd(p,c) (Eq.2).
    
    fd(p,c) = sum of TF-IDF scores of GO terms shared between protein and cluster
    
    Args:
        protein: Protein ID
        cluster: Set of protein IDs in cluster
        protein_go_terms: Dict mapping protein ID to GO terms
        go_tfidf: GOTFIDF instance for TF-IDF scores
        cluster_id: Cluster ID
        
    Returns:
        Functional dependency score
    """
    if protein not in protein_go_terms:
        return 0.0
    
    protein_terms = protein_go_terms[protein]
    if not protein_terms:
        return 0.0
    
    # Sum TF-IDF scores of GO terms that the protein has
    fd_score = 0.0
    for go_term in protein_terms:
        tfidf_score = go_tfidf.get_tfidf(cluster_id, go_term)
        fd_score += tfidf_score
    
    # Normalize by number of protein's GO terms
    if len(protein_terms) > 0:
        fd_score = fd_score / len(protein_terms)
    
    return fd_score


def calculate_membership(protein: str, cluster: Set[str], cluster_id: int,
                         graph: nx.Graph,
                         protein_go_terms: Dict[str, Set[str]],
                         go_tfidf: 'GOTFIDF',
                         permanence_scores: Dict[str, Dict[int, float]],
                         alpha: float = 0.5) -> float:
    """
    Calculate Membership(p,c) (Eq.4).
    
    Membership(p,c) = alpha * Perm(p,c) + (1 - alpha) * fd(p,c)
    
    Args:
        protein: Protein ID
        cluster: Set of protein IDs in cluster
        cluster_id: Cluster ID
        graph: NetworkX graph
        protein_go_terms: Dict mapping protein ID to GO terms
        go_tfidf: GOTFIDF instance
        permanence_scores: Pre-computed permanence scores
        alpha: Weight parameter (0-1) balancing permanence and functional dependency
        
    Returns:
        Membership score
    """
    # Get permanence
    perm = permanence_scores.get(protein, {}).get(cluster_id, 0.0)
    
    # Get functional dependency
    fd = calculate_functional_dependency(protein, cluster, protein_go_terms, 
                                          go_tfidf, cluster_id)
    
    # Calculate membership (Eq.4)
    membership = alpha * perm + (1 - alpha) * fd
    
    return membership


def calculate_intra_extra_links(protein: str, cluster: Set[str], 
                                graph: nx.Graph) -> Tuple[int, int]:
    """
    Calculate intra-cluster and extra-cluster links for a protein.
    
    Args:
        protein: Protein ID
        cluster: Set of protein IDs in cluster
        graph: NetworkX graph
        
    Returns:
        (intra_links, extra_links) tuple
    """
    if protein not in graph:
        return (0, 0)
    
    neighbors = set(graph.neighbors(protein))
    intra_links = len(neighbors & cluster)
    extra_links = len(neighbors - cluster)
    
    return (intra_links, extra_links)


def find_emax_cluster(protein: str, clusters: Dict[int, Set[str]], 
                     graph: nx.Graph) -> int:
    """
    Find the cluster with maximum external connections (E_max).
    
    Args:
        protein: Protein ID
        clusters: Dict mapping cluster_id to set of proteins
        graph: NetworkX graph
        
    Returns:
        Cluster ID with maximum connections
    """
    if protein not in graph:
        return -1
    
    neighbors = set(graph.neighbors(protein))
    max_connections = 0
    emax_cluster_id = -1
    
    for cluster_id, cluster in clusters.items():
        if protein in cluster:
            continue  # Skip the cluster the protein is already in
        
        connections = len(neighbors & cluster)
        if connections > max_connections:
            max_connections = connections
            emax_cluster_id = cluster_id
    
    return emax_cluster_id


def apply_overlap_reassignment(clusters: Dict[int, Set[str]],
                               graph: nx.Graph,
                               protein_go_terms: Dict[str, Set[str]],
                               go_tfidf: 'GOTFIDF',
                               permanence_scores: Dict[str, Dict[int, float]],
                               alpha: float = 0.5,
                               overlap_tau: float = 0.1,
                               transfer_tau: float = 0.0) -> Dict[int, Set[str]]:
    """
    Apply overlapping community reassignment based on membership scores.
    
    Rules:
    1. A protein can be added to additional clusters if membership gain > overlap_tau
    2. A protein can be transferred if Extra-link > Intra-link and transfer_tau threshold met
    
    Args:
        clusters: Initial clusters
        graph: NetworkX graph
        protein_go_terms: Dict mapping protein ID to GO terms
        go_tfidf: GOTFIDF instance
        permanence_scores: Pre-computed permanence scores
        alpha: Membership weight parameter
        overlap_tau: Minimum membership gain to allow overlap
        transfer_tau: Threshold for transfer (Extra-link > Intra-link)
        
    Returns:
        Updated clusters with overlaps
    """
    updated_clusters = {cid: proteins.copy() for cid, proteins in clusters.items()}
    
    # Get all proteins
    all_proteins = set()
    for cluster in clusters.values():
        all_proteins.update(cluster)
    
    logger.info(f"Applying overlap reassignment for {len(all_proteins)} proteins...")
    
    for protein in all_proteins:
        # Find current cluster(s)
        current_clusters = [cid for cid, cluster in updated_clusters.items() 
                           if protein in cluster]
        
        if not current_clusters:
            continue
        
        # Calculate membership in current clusters
        current_memberships = {}
        for cid in current_clusters:
            cluster = updated_clusters[cid]
            memb = calculate_membership(protein, cluster, cid, graph, 
                                       protein_go_terms, go_tfidf, 
                                       permanence_scores, alpha)
            current_memberships[cid] = memb
        
        # Check all other clusters for potential overlap
        for cluster_id, cluster in updated_clusters.items():
            if cluster_id in current_clusters:
                continue
            
            # Calculate membership if added to this cluster
            test_cluster = cluster | {protein}
            memb_if_added = calculate_membership(protein, test_cluster, cluster_id,
                                               graph, protein_go_terms, go_tfidf,
                                               permanence_scores, alpha)
            
            # Check overlap condition: membership gain > threshold
            max_current_memb = max(current_memberships.values()) if current_memberships else 0.0
            membership_gain = memb_if_added - max_current_memb
            
            if membership_gain > overlap_tau:
                # Add protein to this cluster (overlap)
                updated_clusters[cluster_id].add(protein)
                logger.debug(f"Added {protein} to cluster {cluster_id} (gain={membership_gain:.3f})")
        
        # Check transfer condition: Extra-link > Intra-link
        for cid in current_clusters:
            cluster = updated_clusters[cid]
            intra_links, extra_links = calculate_intra_extra_links(protein, cluster, graph)
            
            if extra_links > intra_links:
                # Find cluster with maximum external connections
                emax_cid = find_emax_cluster(protein, updated_clusters, graph)
                
                if emax_cid != -1 and emax_cid != cid:
                    # Check transfer threshold
                    emax_cluster = updated_clusters[emax_cid]
                    emax_intra, emax_extra = calculate_intra_extra_links(protein, emax_cluster, graph)
                    
                    if emax_intra > intra_links:  # Transfer improves intra-links
                        # Transfer protein
                        updated_clusters[cid].discard(protein)
                        updated_clusters[emax_cid].add(protein)
                        logger.debug(f"Transferred {protein} from cluster {cid} to {emax_cid}")
    
    return updated_clusters

