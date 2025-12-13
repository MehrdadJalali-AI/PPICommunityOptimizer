# Troubleshooting Guide

## Common Errors and Solutions

### Error: "GO file not found" in Step 2

**Symptoms:**
```
[Step 2] Loading GO annotations...
ERROR - GO file not found: cache/goa_saccharomyces.gaf.gz
```

**Causes:**
1. Running STRING mode without GOA GAF file
2. Wrong file path specified
3. File not downloaded

**Solutions:**

#### For STRING Mode:
1. **Download GOA GAF file:**
   ```bash
   python download_goa.py
   ```
   Or manually from: https://www.ebi.ac.uk/GOA/downloads

2. **Use Gavin mode instead** (if you have GO.txt):
   ```bash
   python main.py \
     --mode gavin \
     --ppi gavin2006_socioaffinities_rescaled.txt \
     --go-file GO.txt \
     --go-use-symbol \
     --go-taxid 559292 \
     --outdir outputs_gavin/
   ```

3. **Skip STRING mode** (run_all.py will skip automatically):
   ```bash
   python run_all.py
   ```

#### For Gavin Mode:
- Check that `GO.txt` exists in current directory
- Verify file path is correct
- Check file permissions

### Error: "No GO annotations loaded"

**Symptoms:**
```
WARNING - No GO annotations loaded!
```

**Causes:**
1. Wrong `taxid` filter
2. Wrong `use_symbol` setting
3. File format mismatch

**Solutions:**

1. **Check taxid:**
   - For SGD GO: use `--go-taxid 559292` (S. cerevisiae)
   - For GOA: may need different taxid or omit filter

2. **Check use_symbol:**
   - For SGD GAF files: use `--go-use-symbol`
   - For GOA GAF files: don't use `--go-use-symbol`

3. **Verify file format:**
   ```bash
   head -5 GO.txt
   ```
   Should show GAF format with tab-separated columns

### Error: Encoding Issues

**Symptoms:**
```
UnicodeDecodeError: 'utf-8' codec can't decode...
```

**Solutions:**
- The loader tries multiple encodings automatically
- If still failing, check file encoding:
  ```bash
  file -I GO.txt
  ```

### Error: MCL Not Found

**Symptoms:**
```
WARNING - MCL not found. Using NetworkX-based approximation.
```

**Solutions:**
- This is not an error - the pipeline uses Louvain as fallback
- To use real MCL, install:
  ```bash
  # macOS
  brew install mcl
  
  # Linux
  sudo apt-get install mcl
  ```

### Error: Missing Libraries

**Symptoms:**
```
ModuleNotFoundError: No module named 'cdlib'
```

**Solutions:**
```bash
pip install -r requirements.txt
pip install cdlib python-igraph markov-clustering
```

## Quick Diagnostic

Run this to check your setup:

```bash
# Check files
ls -lh GO.txt gavin2006_socioaffinities_rescaled.txt 4932.protein.links.detailed.v11.5.txt

# Test GO loading
python3 -c "from src.go_loader import GOLoader; loader = GOLoader(); result = loader.load_from_gaf('GO.txt', taxid=559292, use_symbol=True); print(f'✓ Loaded {len(result)} proteins')"

# Test Gavin loading
python3 -c "from src.gavin_loader import GavinLoader; loader = GavinLoader(); graph = loader.load('gavin2006_socioaffinities_rescaled.txt'); print(f'✓ Loaded {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges')"
```

## Still Having Issues?

1. Check logs for detailed error messages
2. Verify file formats match expected formats
3. Ensure all required files are present
4. Check Python version (requires 3.7+)

