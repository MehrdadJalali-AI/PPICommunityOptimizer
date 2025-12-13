# Fix: GOA GAF File Download Issue

## Problem

The automatic download scripts are failing because the GOA GAF file URLs have changed or require authentication.

**Error:** `ERROR - GO file not found: cache/goa_saccharomyces.gaf.gz`

## Solution Options

### Option 1: Manual Download (Recommended)

The GOA files are now hosted at EBI. Download manually:

1. **Visit the EBI GOA download page:**
   - Go to: https://www.ebi.ac.uk/GOA/downloads
   - Or direct FTP: https://ftp.ebi.ac.uk/pub/databases/GO/goa/

2. **Find S. cerevisiae file:**
   - Look for `SACCHAROMYCES` directory
   - Download `goa_saccharomyces.gaf.gz` (or latest version)

3. **Save to cache directory:**
   ```bash
   mkdir -p cache
   # Move downloaded file to:
   mv ~/Downloads/goa_saccharomyces.gaf.gz cache/
   ```

### Option 2: Use Alternative GO Source

If GOA is unavailable, you can use the SGD GO file (`GO.txt`) that you already have, but you'll need to adapt it for STRING protein IDs.

### Option 3: Skip STRING Mode

The pipeline will automatically skip STRING mode if the GO file is missing. You can still run:
- ✅ Gavin dataset pipeline (uses `GO.txt`)
- ✅ Method comparison

## Current Status

The scripts (`run_all.py`, `run_all.sh`) have been updated to:
- ✅ Check for GO file before running STRING mode
- ✅ Skip STRING mode gracefully if file is missing
- ✅ Continue with other datasets

## Quick Fix

If you just want to run what's available:

```bash
# This will skip STRING mode and run Gavin + comparison
python run_all.py
```

Or run Gavin only:
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

## Note

The GOA GAF file is only needed for STRING mode. If you're primarily using the Gavin dataset, you don't need it.

