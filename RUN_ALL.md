# Running the Complete Pipeline

This guide explains how to run the entire codebase to generate all results.

## Quick Start

Run everything with a single command:

```bash
./run_all.sh
```

This script will:
1. ✅ Run STRING dataset pipeline → `outputs/`
2. ✅ Run Gavin dataset pipeline → `outputs_gavin/`
3. ✅ Run method comparison → `community_detection_comparison.csv`

## Manual Step-by-Step

If you prefer to run steps individually:

### Step 1: STRING Dataset Pipeline

```bash
python main.py \
  --mode string \
  --taxid 4932 \
  --string-mode download \
  --threshold 700 \
  --go-file cache/goa_saccharomyces.gaf.gz \
  --outdir outputs/ \
  --mcl-inflation 2.0 \
  --alpha 0.5 \
  --overlap-tau 0.1 \
  --transfer-tau 0.0 \
  --lea-population 30 \
  --lea-evaluations 500 \
  --lambda-inter 1.0 \
  --lambda-fragment 0.5 \
  --random-seed 42
```

**Outputs in `outputs/`:**
- `clusters_initial_mcl.csv` - Initial MCL clusters
- `clusters_optimized_lea.csv` - LEA-optimized clusters
- `go_term_importance.csv` - TF-IDF scores for GO terms
- `protein_membership.csv` - Membership scores for all proteins
- `overlap_summary.csv` - Overlapping community assignments
- `evaluation_results.csv` - Evaluation metrics

### Step 2: Gavin Dataset Pipeline

```bash
python main.py \
  --mode gavin \
  --ppi gavin2006_socioaffinities_rescaled.txt \
  --go-file GO.txt \
  --go-use-symbol \
  --go-taxid 559292 \
  --outdir outputs_gavin/ \
  --mcl-inflation 2.0 \
  --alpha 0.5 \
  --overlap-tau 0.1 \
  --transfer-tau 0.0 \
  --lea-population 30 \
  --lea-evaluations 500 \
  --lambda-inter 1.0 \
  --lambda-fragment 0.5 \
  --random-seed 42
```

**Outputs in `outputs_gavin/`:**
- Same CSV files as STRING dataset

### Step 3: Method Comparison

```bash
python compare_methods.py --lea-evaluations 500
```

**Output:**
- `community_detection_comparison.csv` - Comparison of all methods

### Faster Testing (Fewer LEA Evaluations)

For faster runs during testing:

```bash
# Skip LEA optimization entirely
python main.py ... --skip-lea

# Or reduce LEA evaluations
python main.py ... --lea-evaluations 50
python compare_methods.py --lea-evaluations 50
```

## Required Files

### For STRING Mode:
- `4932.protein.links.detailed.v11.5.txt` (or in `cache/`)
- `cache/goa_saccharomyces.gaf.gz` (or similar GOA GAF file)

### For Gavin Mode:
- `gavin2006_socioaffinities_rescaled.txt`
- `GO.txt` (SGD GAF v2.0)

## Expected Runtime

- **STRING Pipeline**: ~5-10 minutes (depends on LEA evaluations)
- **Gavin Pipeline**: ~2-5 minutes
- **Method Comparison**: ~5-15 minutes (depends on LEA evaluations)

**Total**: ~15-30 minutes for full run

## Output Files Summary

### Pipeline Outputs (per dataset):
1. **clusters_initial_mcl.csv** - Initial MCL clustering results
2. **clusters_optimized_lea.csv** - LEA-optimized overlapping communities
3. **go_term_importance.csv** - GO term TF-IDF importance scores
4. **protein_membership.csv** - Detailed membership scores (permanence, fd, membership)
5. **overlap_summary.csv** - Summary of overlapping assignments
6. **evaluation_results.csv** - Evaluation metrics (modularity, conductance, etc.)

### Comparison Output:
- **community_detection_comparison.csv** - Comparison across 8 methods on 2 datasets

## Troubleshooting

### Missing Files
If files are missing, the script will skip that dataset and continue.

### Memory Issues
For large networks, reduce LEA evaluations:
```bash
--lea-evaluations 100
```

### Library Errors
Install missing libraries:
```bash
pip install -r requirements.txt
pip install cdlib python-igraph markov-clustering
```

## Verification

After running, verify outputs:

```bash
# Check STRING outputs
ls -lh outputs/*.csv

# Check Gavin outputs  
ls -lh outputs_gavin/*.csv

# Check comparison
ls -lh community_detection_comparison.csv
```

## Next Steps

After generating results:
1. Analyze CSV files in your preferred tool (Python, R, Excel)
2. Visualize results using the evaluation metrics
3. Compare methods using `community_detection_comparison.csv`
4. Update manuscript with new results

