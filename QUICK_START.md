# Quick Start Guide

## Run Everything (Easiest Method)

### Option 1: Bash Script (Linux/macOS)
```bash
./run_all.sh
```

### Option 2: Python Script (Cross-platform)
```bash
python run_all.py
```

Both scripts will:
1. ✅ Run STRING dataset pipeline → `outputs/`
2. ✅ Run Gavin dataset pipeline → `outputs_gavin/`
3. ✅ Run method comparison → `community_detection_comparison.csv`

## Individual Commands

### STRING Dataset Only
```bash
python main.py \
  --mode string \
  --taxid 4932 \
  --string-mode download \
  --threshold 700 \
  --go-file cache/goa_saccharomyces.gaf.gz \
  --outdir outputs/ \
  --random-seed 42
```

### Gavin Dataset Only
```bash
python main.py \
  --mode gavin \
  --ppi gavin2006_socioaffinities_rescaled.txt \
  --go-file GO.txt \
  --go-use-symbol \
  --go-taxid 559292 \
  --outdir outputs_gavin/ \
  --random-seed 42
```

### Method Comparison Only
```bash
python compare_methods.py
```

## Faster Testing (Fewer Evaluations)

```bash
# Reduce LEA evaluations for faster runs
python main.py ... --lea-evaluations 50
python compare_methods.py --lea-evaluations 50

# Or skip LEA entirely
python main.py ... --skip-lea
```

## Expected Outputs

### Per Dataset (`outputs/` or `outputs_gavin/`):
- `clusters_initial_mcl.csv` - Initial MCL clusters
- `clusters_optimized_lea.csv` - LEA-optimized clusters
- `go_term_importance.csv` - GO term TF-IDF scores
- `protein_membership.csv` - Membership scores
- `overlap_summary.csv` - Overlapping assignments
- `evaluation_results.csv` - Evaluation metrics

### Comparison:
- `community_detection_comparison.csv` - All methods comparison

## Runtime Estimates

- **Full pipeline (both datasets)**: ~15-30 minutes
- **Single dataset**: ~5-10 minutes
- **Method comparison**: ~5-15 minutes
- **Fast mode (50 evaluations)**: ~2-5 minutes per dataset

## Troubleshooting

**Missing files?** Scripts will skip missing datasets and continue.

**Memory issues?** Reduce `--lea-evaluations` to 100 or 50.

**Library errors?** Install dependencies:
```bash
pip install -r requirements.txt
pip install cdlib python-igraph markov-clustering
```

## More Details

See [RUN_ALL.md](RUN_ALL.md) for detailed instructions.

