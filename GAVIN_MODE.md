# Gavin PPI + SGD GO Mode

This document describes how to run the community detection pipeline using Gavin PPI network and SGD GO annotations.

## Data Sources

- **Gavin PPI**: `gavin2006_socioaffinities_rescaled.txt` - Weighted socio-affinity protein-protein interactions
- **SGD GO**: `GO.txt` - Gene Ontology annotations from Saccharomyces Genome Database (GAF v2.0 format)

## Quick Start

```bash
python main.py \
  --mode gavin \
  --ppi gavin2006_socioaffinities_rescaled.txt \
  --go-file GO.txt \
  --go-use-symbol \
  --go-taxid 559292 \
  --alpha 0.5 \
  --overlap-tau 0.1 \
  --outdir outputs_gavin/
```

Or use the provided script:

```bash
./run_gavin.sh
```

## File Formats

### Gavin PPI Format

Tab-separated file with header:
```
s	d	description
YKL144C	YPR110C	0.397689
YKL144C	YOR210W	0.253468
...
```

- Column 1: Source protein (yeast ORF name)
- Column 2: Target protein (yeast ORF name)
- Column 3: Socio-affinity weight (0-1, normalized)

### SGD GO Format (GAF v2.0)

Tab-separated GAF file:
```
SGD	S000001503	SPT23		GO:0003674	...
SGD	S000000735	GCN4		GO:1990139	...
```

- Column 2: DB_Object_Symbol (yeast ORF name, e.g., YDL159W)
- Column 5: GO_ID (e.g., GO:0003674)
- Use `--go-use-symbol` to extract protein IDs from DB_Object_Symbol column

## Key Differences from STRING Mode

1. **PPI Source**: Gavin socio-affinity vs STRING combined scores
2. **Protein IDs**: Yeast ORF names (YDL159W) vs STRING/UniProt IDs
3. **GO Format**: SGD GAF v2.0 vs GOA GAF format
4. **Weight Normalization**: Gavin weights are already normalized [0,1]

## Output Files

Same structure as STRING mode:
- `clusters_initial_mcl.csv`
- `go_term_importance.csv`
- `protein_membership.csv`
- `clusters_optimized_lea.csv`
- `overlap_summary.csv`
- `evaluation_results.csv`

## Parameters

- `--mode gavin`: Use Gavin PPI loader
- `--ppi`: Path to Gavin PPI file
- `--go-file`: Path to SGD GO file
- `--go-use-symbol`: Use DB_Object_Symbol for protein IDs (required for SGD)
- `--go-taxid 559292`: Filter GO annotations for yeast (S. cerevisiae)

## Network Statistics

- **Nodes**: ~1,860 proteins
- **Edges**: ~7,601 interactions
- **GO Annotations**: ~6,439 proteins with GO terms

## Notes

- Protein IDs must match between PPI and GO files (both use yeast ORF names)
- The pipeline automatically normalizes Gavin weights if needed
- MCL clustering works on weighted adjacency matrix
- All equations (permanence, fd, membership) work identically to STRING mode

