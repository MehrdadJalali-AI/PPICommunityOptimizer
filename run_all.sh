#!/bin/bash
# Master script to run entire pipeline and generate all results
# This script runs:
# 1. STRING dataset pipeline
# 2. Gavin dataset pipeline  
# 3. Method comparison
# 4. Generates all CSV outputs

set -e  # Exit on error

echo "=========================================="
echo "Community Detection Pipeline - Full Run"
echo "=========================================="
echo ""

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if required files exist
echo -e "${BLUE}Checking required files...${NC}"

if [ ! -f "4932.protein.links.detailed.v11.5.txt" ] && [ ! -f "cache/4932.protein.links.detailed.v11.5.txt" ] && [ ! -f "cache/4932.protein.links.detailed.v11.5.txt.gz" ]; then
    echo -e "${YELLOW}Warning: STRING PPI file not found. STRING mode will be skipped.${NC}"
    echo "  Expected: 4932.protein.links.detailed.v11.5.txt (or in cache/)"
    SKIP_STRING=true
else
    SKIP_STRING=false
fi

if [ ! -f "gavin2006_socioaffinities_rescaled.txt" ]; then
    echo -e "${YELLOW}Warning: Gavin PPI file not found. Gavin mode will be skipped.${NC}"
    SKIP_GAVIN=true
else
    SKIP_GAVIN=false
fi

if [ ! -f "GO.txt" ]; then
    echo -e "${YELLOW}Warning: GO.txt not found. Gavin mode will be skipped.${NC}"
    SKIP_GAVIN=true
else
    SKIP_GAVIN=false
fi

echo ""

# Step 1: STRING Dataset Pipeline
if [ "$SKIP_STRING" = false ]; then
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Step 1: Running STRING Dataset Pipeline${NC}"
    echo -e "${GREEN}========================================${NC}"
    
    # Find GO file
    GO_FILE=""
    if [ -f "cache/goa_saccharomyces.gaf.gz" ]; then
        GO_FILE="cache/goa_saccharomyces.gaf.gz"
    elif [ -f "cache/goa_saccharomyces.gaf" ]; then
        GO_FILE="cache/goa_saccharomyces.gaf"
    else
        echo -e "${YELLOW}Warning: GOA GAF file not found for STRING mode.${NC}"
        echo "  To download: wget http://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz -P cache/"
        echo -e "${YELLOW}Skipping STRING dataset (GO file required).${NC}"
        SKIP_STRING=true
    fi
    
    python main.py \
        --mode string \
        --taxid 4932 \
        --string-mode download \
        --threshold 700 \
        --go-file "$GO_FILE" \
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
    
    echo -e "${GREEN}✓ STRING pipeline completed${NC}"
    echo ""
else
    echo -e "${YELLOW}Skipping STRING dataset (files not found)${NC}"
    echo ""
fi

# Step 2: Gavin Dataset Pipeline
if [ "$SKIP_GAVIN" = false ]; then
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}Step 2: Running Gavin Dataset Pipeline${NC}"
    echo -e "${GREEN}========================================${NC}"
    
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
    
    echo -e "${GREEN}✓ Gavin pipeline completed${NC}"
    echo ""
else
    echo -e "${YELLOW}Skipping Gavin dataset (files not found)${NC}"
    echo ""
fi

# Step 3: Method Comparison
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Step 3: Running Method Comparison${NC}"
echo -e "${GREEN}========================================${NC}"

python compare_methods.py --lea-evaluations 500

echo -e "${GREEN}✓ Method comparison completed${NC}"
echo ""

# Summary
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Summary of Generated Files${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

if [ "$SKIP_STRING" = false ]; then
    echo "STRING Dataset Outputs (outputs/):"
    ls -lh outputs/*.csv 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
    echo ""
fi

if [ "$SKIP_GAVIN" = false ]; then
    echo "Gavin Dataset Outputs (outputs_gavin/):"
    ls -lh outputs_gavin/*.csv 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
    echo ""
fi

echo "Comparison Results:"
ls -lh community_detection_comparison.csv 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}All pipelines completed successfully!${NC}"
echo -e "${GREEN}========================================${NC}"

