# Fix: STRING Dataset Jaccard Calculation

## Problem

Jaccard similarity is **empty for all STRING methods** because:
1. The GO annotation file (`cache/goa_saccharomyces.gaf.gz`) is **not downloaded**
2. Without GO annotations, `protein_go_terms` is empty
3. `calculate_go_jaccard_similarity()` returns 0.0 when no GO terms are available

## Root Cause

The `load_string_dataset()` function in `compare_methods.py` was looking for:
- `cache/goa_4932.gaf.gz` (taxid-based naming)

But `download_goa.py` downloads:
- `cache/goa_saccharomyces.gaf.gz` (organism-based naming)

**This mismatch has been fixed** - the code now checks for both filenames.

## Solution

### Step 1: Download GO Annotations

Run the download script:

```bash
python download_goa.py
```

This will download `cache/goa_saccharomyces.gaf.gz` (~40-50 MB).

### Step 2: Regenerate Results

After downloading, regenerate the CSV files:

```bash
python generate_updated_results.py
```

This will:
- Load GO annotations for STRING dataset
- Calculate Jaccard similarity for all methods
- Generate updated CSV files with `mean_go_jaccard` populated

## Verification

After downloading and regenerating, you should see:

1. **GO terms loaded**: 
   ```python
   from compare_methods import load_string_dataset
   graph_str, lea_data_str = load_string_dataset()
   print(f"GO terms: {len(lea_data_str.get('protein_go_terms', {}))}")
   # Should show: GO terms: ~6000+ (not 0)
   ```

2. **Jaccard values in CSV**:
   - `results_string_updated.csv` should have `mean_go_jaccard` values (not empty)
   - Values should be in range [0, 1]

## Code Changes Made

✅ **Fixed filename mismatch** in `compare_methods.py`:
- Now checks for multiple possible filenames:
  - `cache/goa_saccharomyces.gaf.gz` (standard)
  - `cache/goa_4932.gaf.gz` (alternative)
  - Uncompressed versions (.gaf)
- Added logging to show which file was found
- Added warning if no GO file is found

## Expected Results

After fixing:

| Dataset | Method | mean_go_jaccard |
|---------|--------|-----------------|
| STRING  | Louvain | ~0.15-0.25 |
| STRING  | Leiden  | ~0.15-0.25 |
| STRING  | MCL     | ~0.15-0.25 |
| STRING  | LEA     | ~0.15-0.25 |
| ...     | ...     | ... |

(Exact values depend on the communities detected)

## Notes

- **Gavin dataset**: Already has GO annotations loaded (from `GO.txt`), so Jaccard works there
- **STRING dataset**: Requires download of GOA GAF file
- The Jaccard calculation itself is working correctly - it just needs GO data

