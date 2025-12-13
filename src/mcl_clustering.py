"""
MCL (Markov Cluster Algorithm) clustering for initial community detection.
"""

import logging
import subprocess
import tempfile
import os
import networkx as nx
from typing import Dict, List, Set

logger = logging.getLogger(__name__)


class MCLClustering:
    """
    MCL clustering wrapper.
    Requires MCL to be installed: https://micans.org/mcl/
    """
    
    def __init__(self, inflation: float = 2.0, min_cluster_size: int = 10):
        """
        Initialize MCL clustering.
        
        Args:
            inflation: MCL inflation parameter (higher = more clusters)
            min_cluster_size: Minimum cluster size to include (filters small clusters)
                            Default: 10 (matches biological interpretability standards)
        """
        self.inflation = inflation
        self.min_cluster_size = min_cluster_size
    
    def cluster(self, graph: nx.Graph) -> Dict[int, Set[str]]:
        """
        Run MCL clustering on graph.
        
        Args:
            graph: NetworkX graph
            
        Returns:
            clusters: Dict mapping cluster_id to set of protein IDs
        """
        logger.info(f"Running MCL clustering (inflation={self.inflation})...")
        
        # Check if MCL is available
        try:
            subprocess.run(['mcl', '--version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("MCL not found. Using NetworkX-based approximation.")
            return self._fallback_clustering(graph)
        
        # Write graph to temporary file in MCL format
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp_in:
            # Write edges with weights
            for u, v, data in graph.edges(data=True):
                weight = data.get('weight', 1.0)
                tmp_in.write(f"{u}\t{v}\t{weight}\n")
            tmp_in_path = tmp_in.name
        
        tmp_out_path = tmp_in_path + '.out'
        
        try:
            # Run MCL
            cmd = ['mcl', tmp_in_path, '--abc', '-I', str(self.inflation), '-o', tmp_out_path]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse MCL output and filter small clusters
            clusters = {}
            cluster_id = 0
            filtered_count = 0
            
            with open(tmp_out_path, 'r') as f:
                for line in f:
                    proteins = line.strip().split('\t')
                    if proteins and proteins[0]:
                        cluster_set = set(proteins)
                        # Filter clusters smaller than min_cluster_size
                        if len(cluster_set) >= self.min_cluster_size:
                            clusters[cluster_id] = cluster_set
                            cluster_id += 1
                        else:
                            filtered_count += 1
            
            logger.info(f"MCL found {len(clusters)} clusters (filtered {filtered_count} clusters < {self.min_cluster_size} proteins)")
            return clusters
            
        except subprocess.CalledProcessError as e:
            logger.error(f"MCL failed: {e.stderr}")
            return self._fallback_clustering(graph)
        finally:
            # Cleanup
            if os.path.exists(tmp_in_path):
                os.remove(tmp_in_path)
            if os.path.exists(tmp_out_path):
                os.remove(tmp_out_path)
    
    def _fallback_clustering(self, graph: nx.Graph) -> Dict[int, Set[str]]:
        """
        Fallback to NetworkX community detection if MCL is not available.
        Uses Louvain algorithm as approximation.
        """
        logger.info("Using Louvain algorithm as fallback...")
        try:
            import community.community_louvain as community_louvain
            partition = community_louvain.best_partition(graph)
            
            clusters = {}
            filtered_count = 0
            
            for node, cluster_id in partition.items():
                if cluster_id not in clusters:
                    clusters[cluster_id] = set()
                clusters[cluster_id].add(node)
            
            # Filter small clusters
            filtered_clusters = {}
            cluster_id_new = 0
            for old_id, cluster_set in clusters.items():
                if len(cluster_set) >= self.min_cluster_size:
                    filtered_clusters[cluster_id_new] = cluster_set
                    cluster_id_new += 1
                else:
                    filtered_count += 1
            
            logger.info(f"Louvain found {len(filtered_clusters)} clusters (filtered {filtered_count} clusters < {self.min_cluster_size} proteins)")
            return filtered_clusters
        except ImportError:
            logger.warning("python-louvain not available. Using simple connected components.")
            clusters = {}
            for i, component in enumerate(nx.connected_components(graph)):
                clusters[i] = component
            logger.info(f"Found {len(clusters)} connected components")
            return clusters

