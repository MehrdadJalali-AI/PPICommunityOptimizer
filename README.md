# Community Detection Pipeline with LEA Optimization

A complete Python pipeline for protein-protein interaction (PPI) network community detection using the Lotus Effect Algorithm (LEA) for optimization.

## Overview

This pipeline implements:
1. **STRING PPI Network Loading**: Download or API-based loading of protein-protein interactions
2. **GO Annotation Loading**: Gene Ontology annotations from GOA GAF files
3. **MCL Clustering**: Initial community detection using Markov Cluster Algorithm
4. **GO TF-IDF**: Term importance calculation for GO terms in clusters (Eq.3)
5. **Permanence Calculation**: Node permanence in communities (Eq.1)
6. **Functional Dependency**: GO-based functional similarity (Eq.2)
7. **Membership Calculation**: Combined permanence and functional dependency (Eq.4)
8. **Overlapping Communities**: Support for proteins in multiple clusters
9. **LEA Optimization**: Optimization of membership parameters using Lotus Effect Algorithm
10. **Evaluation**: Comprehensive metrics including intra/inter density, conductance, modularity, and biological coherence

## Installation

### Prerequisites

- Python 3.7+
- MCL (Markov Cluster Algorithm) - optional but recommended

### Install MCL

**macOS:**
```bash
brew install mcl
```

**Linux:**
```bash
sudo apt-get install mcl
# or
sudo yum install mcl
```

**From source:**
```bash
wget https://micans.org/mcl/src/mcl-14-137.tar.gz
tar -xzf mcl-14-137.tar.gz
cd mcl-14-137
./configure
make
sudo make install
```

### Install Python Dependencies

```bash
pip install -r requirements.txt
```

## Data Preparation

### STRING Database Files

Download STRING database files for your organism:

1. **Protein Links** (detailed):
   ```
   https://stringdb-static.org/download/protein.links.detailed.v11.5/{taxid}.protein.links.detailed.v11.5.txt.gz
   ```

2. **Protein Aliases**:
   ```
   https://stringdb-static.org/download/protein.aliases.v11.5/{taxid}.protein.aliases.v11.5.txt.gz
   ```

Place these files in the `cache/` directory (or specify with `--cache-dir`).

Example for S. cerevisiae (taxid=4932):
```bash
mkdir -p cache
cd cache
wget https://stringdb-static.org/download/protein.links.detailed.v11.5/4932.protein.links.detailed.v11.5.txt.gz
wget https://stringdb-static.org/download/protein.aliases.v11.5/4932.protein.aliases.v11.5.txt.gz
```

### GO Annotations (GOA GAF)

Download GOA GAF file for your organism:

```bash
cd cache
wget http://geneontology.org/gene-associations/goa_{organism}.gaf.gz
```

For S. cerevisiae:
```bash
wget http://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz
```

## Usage

### Basic Usage

```bash
python main.py --taxid 4932 --string-mode download --threshold 700 --go-source goa --outdir outputs/
```

### Command-Line Arguments

- `--taxid`: NCBI taxonomy ID (required)
  - Example: `4932` for S. cerevisiae, `9606` for H. sapiens
  
- `--string-mode`: STRING data loading mode
  - `download`: Parse downloaded files (default, recommended)
  - `api`: Use STRING REST API (limited)
  
- `--threshold`: STRING combined score threshold (0-1000)
  - Default: `700` (high confidence)
  - Lower values (e.g., `400`) include more interactions
  
- `--go-source`: GO annotation source
  - Currently only `goa` supported
  
- `--go-file`: Path to GOA GAF file
  - If not provided, looks for `cache/goa_{taxid}.gaf.gz`
  
- `--outdir`: Output directory for results
  - Default: `outputs/`
  
- `--cache-dir`: Cache directory for downloaded files
  - Default: `cache/`
  
- `--mcl-inflation`: MCL inflation parameter
  - Default: `2.0`
  - Higher values create more, smaller clusters
  
- `--lea-population`: LEA population size
  - Default: `30`
  
- `--lea-evaluations`: Maximum LEA function evaluations
  - Default: `500`
  - Increase for better optimization (slower)
  
- `--lambda-inter`: Weight for inter-cluster penalty
  - Default: `1.0`
  
- `--lambda-fragment`: Weight for fragmentation penalty
  - Default: `0.5`
  
- `--random-seed`: Random seed for reproducibility
  - Example: `42`
  
- `--skip-lea`: Skip LEA optimization (use initial MCL clusters)
  
- `--gold-standard`: Path to gold standard clusters CSV (optional)
  - Format: `cluster_id,protein_id`

### Example Commands

**S. cerevisiae with default settings:**
```bash
python main.py --taxid 4932 --string-mode download --threshold 700 --go-source goa --outdir outputs/
```

**H. sapiens with custom parameters:**
```bash
python main.py --taxid 9606 --string-mode download --threshold 600 --mcl-inflation 2.5 --lea-evaluations 1000 --random-seed 42 --outdir outputs_human/
```

**Quick test without LEA optimization:**
```bash
python main.py --taxid 4932 --string-mode download --threshold 700 --skip-lea --outdir outputs_quick/
```

## Output Files

All outputs are saved as CSV files in the specified output directory:

1. **`clusters_initial_mcl.csv`**
   - Initial MCL clusters
   - Columns: `cluster_id`, `protein_id`

2. **`go_term_importance.csv`**
   - TF-IDF scores for GO terms in clusters
   - Columns: `cluster_id`, `go_term`, `tfidf_score`

3. **`protein_membership.csv`**
   - Detailed membership metrics for each protein-cluster pair
   - Columns: `protein_id`, `cluster_id`, `permanence`, `fd`, `membership`, `intra`, `extra`, `emax`

4. **`clusters_optimized_lea.csv`**
   - Optimized clusters after LEA
   - Columns: `cluster_id`, `protein_id`, `membership_score`

5. **`overlap_summary.csv`**
   - Summary of overlapping communities
   - Columns: `protein_id`, `num_clusters`, `clusters_json`

6. **`evaluation_results.csv`**
   - Evaluation metrics
   - Columns: `intra_density_mean`, `inter_density_mean`, `conductance_mean`, `overlapping_modularity`, `mean_fd_per_cluster`, `num_clusters`, `num_singletons`, etc.

## Project Structure

```
CommunityDetection/
├── main.py                 # Main CLI entry point
├── requirements.txt        # Python dependencies
├── README.md              # This file
├── src/
│   ├── __init__.py
│   ├── string_loader.py   # STRING PPI network loading
│   ├── go_loader.py       # GO annotation loading
│   ├── mcl_clustering.py  # MCL clustering
│   ├── go_tfidf.py        # GO TF-IDF calculation (Eq.3)
│   ├── permanence.py      # Permanence calculation (Eq.1)
│   ├── membership_overlap.py  # Membership & overlap (Eq.2, Eq.4)
│   ├── evaluation.py     # Evaluation metrics
│   ├── outputs.py         # Output file generation
│   └── lea/
│       ├── __init__.py
│       ├── lotus_effect_algorithm.py  # LEA core algorithm
│       ├── fitness_membership.py      # Fitness function
│       └── optimize.py                # Optimization wrapper
├── cache/                 # Cache directory for downloaded files
├── outputs/              # Output directory for results
└── tests/                # Test files (optional)
```

## Caching Behavior

- Downloaded STRING files are cached in `cache/` directory
- The pipeline checks for existing files before downloading
- GO annotations are loaded from GAF files (not cached separately)

## Reproducibility

Use `--random-seed` to ensure reproducible results:

```bash
python main.py --taxid 4932 --random-seed 42 --outdir outputs_reproducible/
```

## Performance Notes

- **Large networks**: For networks with >10,000 proteins, consider:
  - Increasing `--threshold` to reduce network size
  - Using `--skip-lea` for faster initial exploration
  - Reducing `--lea-evaluations` for quicker optimization

- **Memory usage**: Large networks may require significant RAM
  - STRING files can be several GB uncompressed
  - NetworkX graphs are memory-efficient but scale with edges

## Troubleshooting

### MCL not found
If MCL is not installed, the pipeline falls back to Louvain algorithm (requires `python-louvain` package).

### STRING files not found
Ensure files are downloaded and placed in the cache directory with correct naming:
- `{taxid}.protein.links.detailed.v11.5.txt.gz`
- `{taxid}.protein.aliases.v11.5.txt.gz`

### GO annotations missing
If GO annotations are not available, the pipeline will still run but:
- Functional dependency (fd) will be 0
- Membership will rely only on permanence
- GO coherence metrics will be unavailable

## Citation

If you use this code, please cite:
- The Lotus Effect Algorithm (LEA) from MultiModelLEA.ipynb
- STRING database: https://string-db.org/
- Gene Ontology: http://geneontology.org/

## License

[Specify your license here]

## Contact

[Your contact information]

