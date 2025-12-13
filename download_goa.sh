#!/bin/bash
# Download GOA GAF file for S. cerevisiae

echo "Downloading GOA GAF file for S. cerevisiae..."
echo ""

mkdir -p cache
cd cache

# Try multiple URLs
URLS=(
    "https://ftp.ebi.ac.uk/pub/databases/GO/goa/SACCHAROMYCES/goa_saccharomyces.gaf.gz"
    "http://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz"
    "https://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz"
)

DOWNLOADED=false
for URL in "${URLS[@]}"; do
    echo "Trying: $URL"
    if wget "$URL" -O goa_saccharomyces.gaf.gz; then
        DOWNLOADED=true
        break
    else
        echo "  Failed, trying next URL..."
    fi
done

if [ "$DOWNLOADED" = true ]; then
    echo ""
    echo "✓ Successfully downloaded: cache/goa_saccharomyces.gaf.gz"
    echo ""
    ls -lh goa_saccharomyces.gaf.gz
    echo ""
    echo "You can now run the STRING pipeline:"
    echo "  python run_all.py"
else
    echo ""
    echo "✗ Download failed. Please check your internet connection."
    echo ""
    echo "Alternative: Download manually from:"
    echo "  http://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz"
    echo "  Save to: cache/goa_saccharomyces.gaf.gz"
fi

cd ..

