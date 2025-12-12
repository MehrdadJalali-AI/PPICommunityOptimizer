"""
Sanity test: Run pipeline on a tiny toy graph.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import networkx as nx
from src.mcl_clustering import MCLClustering
from src.go_tfidf import GOTFIDF
from src.permanence import calculate_permanence_all_proteins
from src.evaluation import evaluate_clusters

def test_toy_graph():
    """Test pipeline components on a small toy graph."""
    print("Creating toy PPI graph...")
    
    # Create a small graph with 3 communities
    graph = nx.Graph()
    
    # Community 1: proteins A1-A5
    for i in range(1, 6):
        graph.add_edge(f'A{i}', f'A{i+1}' if i < 5 else 'A1')
    
    # Community 2: proteins B1-B5
    for i in range(1, 6):
        graph.add_edge(f'B{i}', f'B{i+1}' if i < 5 else 'B1')
    
    # Community 3: proteins C1-C5
    for i in range(1, 6):
        graph.add_edge(f'C{i}', f'C{i+1}' if i < 5 else 'C1')
    
    # Add some inter-community edges
    graph.add_edge('A1', 'B1')
    graph.add_edge('B1', 'C1')
    
    print(f"Graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
    
    # Test MCL clustering
    print("\nTesting MCL clustering...")
    mcl = MCLClustering(inflation=2.0)
    clusters = mcl.cluster(graph)
    print(f"Found {len(clusters)} clusters")
    for cid, cluster in clusters.items():
        print(f"  Cluster {cid}: {len(cluster)} proteins")
    
    # Create dummy GO terms
    print("\nTesting GO TF-IDF...")
    protein_go_terms = {}
    for node in graph.nodes():
        # Assign some dummy GO terms
        protein_go_terms[node] = {f'GO:000{i:04d}' for i in range(1, 4)}
    
    go_tfidf = GOTFIDF(clusters, protein_go_terms)
    print("GO TF-IDF computed successfully")
    
    # Test permanence
    print("\nTesting permanence calculation...")
    permanence_scores = calculate_permanence_all_proteins(clusters, graph)
    print(f"Computed permanence for {len(permanence_scores)} proteins")
    
    # Test evaluation
    print("\nTesting evaluation...")
    eval_df = evaluate_clusters(clusters, graph, protein_go_terms, go_tfidf)
    print("Evaluation metrics:")
    print(eval_df.to_string())
    
    print("\n✓ All tests passed!")


if __name__ == '__main__':
    test_toy_graph()

