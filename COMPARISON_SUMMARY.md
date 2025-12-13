# Community Detection Comparison Summary

## Overview

This comparison evaluates the LEA-based overlapping community detection method against state-of-the-art algorithms on PPI networks.

## Methods Compared

### Non-Overlapping Methods:
1. **Louvain** - Fast modularity optimization
2. **Leiden** - Improved Louvain with guarantees
3. **Infomap** - Information-theoretic approach
4. **MCL** - Markov Cluster Algorithm
5. **OSLOM** - Order Statistics Local Optimization Method (approximated via label propagation)

### Overlapping Methods:
6. **Link Communities** - Edge-based overlapping communities
7. **BigCLAM** - Overlapping community detection via NMF
8. **LEA_Overlapping** - Our method with permanence + functional dependency optimization

## Datasets

1. **STRING** - S. cerevisiae PPI network (taxid 4932)
   - Threshold: 700 (combined_score)
   - Nodes: ~6000
   - Edges: ~120,000

2. **Gavin** - Yeast socio-affinity PPI network
   - Normalized socio-affinity scores
   - Nodes: ~2000
   - Edges: ~7000

## Evaluation Metrics

- **Modularity**: Quality of community structure
- **Conductance**: Average boundary quality
- **NMI**: Normalized Mutual Information (if ground truth available)
- **Overlapping NMI**: Overlap-aware NMI
- **Runtime**: Execution time in seconds

## Usage

### Full Comparison (may take hours):
```bash
python compare_methods.py
```

### Fast Test (skip LEA, fewer evaluations):
```bash
python compare_methods.py --skip-lea --lea-evaluations 50
```

### Custom LEA Evaluations:
```bash
python compare_methods.py --lea-evaluations 200
```

## Output

Results are saved to `community_detection_comparison.csv` with columns:
- dataset, method, overlapping
- num_nodes, num_edges, num_communities
- avg_community_size, modularity, nmi, overlapping_nmi
- conductance, runtime_sec, parameters

## Reproducibility

- Fixed random seed: 42
- All parameters logged in JSON format
- Same preprocessing for all methods

