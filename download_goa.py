#!/usr/bin/env python3
"""
Download GOA GAF file for S. cerevisiae.
"""

import os
import sys
import urllib.request
from pathlib import Path

def download_goa():
    """Download GOA GAF file."""
    # Try multiple possible URLs
    urls = [
        "https://ftp.ebi.ac.uk/pub/databases/GO/goa/SACCHAROMYCES/goa_saccharomyces.gaf.gz",
        "http://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz",
        "https://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz"
    ]
    cache_dir = Path("cache")
    output_file = cache_dir / "goa_saccharomyces.gaf.gz"
    
    print("Downloading GOA GAF file for S. cerevisiae...")
    print(f"Destination: {output_file}")
    print()
    
    # Create cache directory
    cache_dir.mkdir(exist_ok=True)
    
    # Check if file already exists
    if output_file.exists():
        print(f"⚠ File already exists: {output_file}")
        response = input("Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print("Skipping download.")
            return True
    
    # Try each URL until one works
    url_used = None
    for url in urls:
        try:
            print(f"Trying: {url}")
            # Download with progress
            def show_progress(block_num, block_size, total_size):
                downloaded = block_num * block_size
                percent = min(100, (downloaded * 100) / total_size) if total_size > 0 else 0
                print(f"\rDownloading... {percent:.1f}%", end='', flush=True)
            
            print("Starting download...")
            urllib.request.urlretrieve(url, output_file, show_progress)
            url_used = url
            break  # Success, exit loop
        except urllib.error.HTTPError as e:
            if e.code == 403 or e.code == 404:
                print(f"  ✗ Failed (HTTP {e.code}), trying next URL...")
                continue
            else:
                raise
        except Exception as e:
            print(f"  ✗ Error: {e}, trying next URL...")
            continue
    if url_used is None:
        # All URLs failed
        print("\n✗ All download URLs failed.")
        print("\nPlease download manually from one of these URLs:")
        for url in urls:
            print(f"  {url}")
        print(f"\nSave to: {output_file}")
        return False
    
    try:
        print("\n")
        
        # Check file size
        file_size = output_file.stat().st_size
        print(f"✓ Successfully downloaded: {output_file}")
        print(f"  Size: {file_size:,} bytes ({file_size / 1024 / 1024:.2f} MB)")
        print()
        print("You can now run the STRING pipeline:")
        print("  python run_all.py")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Error after download: {e}")
        if output_file.exists():
            output_file.unlink()  # Remove partial file
        return False

if __name__ == "__main__":
    success = download_goa()
    sys.exit(0 if success else 1)

