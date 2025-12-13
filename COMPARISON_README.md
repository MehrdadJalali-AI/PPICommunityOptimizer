# Community Detection Method Comparison

This script compares the LEA-based overlapping community detection method against state-of-the-art algorithms on PPI networks.

## Installation

Install required libraries:

```bash
pip install -r requirements.txt

# Additional libraries for comparison:
pip install cdlib python-igraph markov-clustering
```

**Note**: Some algorithms may not be available if libraries are missing. The script will use fallbacks where possible.

## Usage

Run the comparison:

```bash
python compare_methods.py
```

This will:
1. Load STRING dataset (if available)
2. Load Gavin dataset
3. Run all methods on both datasets
4. Generate `community_detection_comparison.csv`

## Methods Compared

### Non-Overlapping Methods:
- **Louvain**: Fast modularity optimization
- **Leiden**: Improved version of Louvain
- **Infomap**: Information-theoretic approach
- **MCL**: Markov Cluster Algorithm
- **OSLOM**: Order Statistics Local Optimization Method (approximated)

### Overlapping Methods:
- **Link Communities**: Edge-based overlapping communities
- **BigCLAM**: Overlapping community detection via NMF
- **LEA_Overlapping**: Our LEA-based method with permanence + functional dependency

## Output Format

The CSV file `community_detection_comparison.csv` contains:

- `dataset`: Dataset name (STRING or Gavin)
- `method`: Algorithm name
- `overlapping`: True/False
- `num_nodes`: Number of nodes in graph
- `num_edges`: Number of edges in graph
- `num_communities`: Number of detected communities
- `avg_community_size`: Average community size
- `modularity`: Modularity score
- `nmi`: Normalized Mutual Information (if ground truth available)
- `overlapping_nmi`: Overlapping NMI (if ground truth available)
- `conductance`: Average conductance
- `runtime_sec`: Runtime in seconds
- `parameters`: JSON string with method parameters

## Reproducibility

- Random seed: 42 (fixed)
- All methods use same graph preprocessing
- Parameters logged in JSON format

## Notes

- OSLOM requires external binary; script uses label propagation as approximation
- Some methods may fail on large networks; errors are logged but don't stop execution
- LEA method requires GO annotations for full functionality

