#!/usr/bin/env python3
"""
Master script to run entire pipeline and generate all results.
Cross-platform Python version.
"""

import os
import sys
import subprocess
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"\n{'='*60}")
    print(f"{description}")
    print(f"{'='*60}")
    print(f"Running: {' '.join(cmd)}")
    print()
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        print(f"\n✓ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n✗ {description} failed with error code {e.returncode}")
        return False
    except FileNotFoundError:
        print(f"\n✗ Command not found. Make sure Python is in PATH.")
        return False

def check_file(filepath, description):
    """Check if a file exists."""
    if os.path.exists(filepath):
        print(f"✓ Found: {filepath}")
        return True
    else:
        print(f"✗ Missing: {filepath} ({description})")
        return False

def main():
    """Run all pipelines."""
    print("="*60)
    print("Community Detection Pipeline - Full Run")
    print("="*60)
    print()
    
    # Check required files
    print("Checking required files...")
    print("-"*60)
    
    string_ppi_exists = (
        check_file("4932.protein.links.detailed.v11.5.txt", "STRING PPI file") or
        check_file("cache/4932.protein.links.detailed.v11.5.txt", "STRING PPI file (cache)") or
        check_file("cache/4932.protein.links.detailed.v11.5.txt.gz", "STRING PPI file (cache, gz)")
    )
    
    gavin_ppi_exists = check_file("gavin2006_socioaffinities_rescaled.txt", "Gavin PPI file")
    go_file_exists = check_file("GO.txt", "GO annotations (Gavin)")
    
    # Find GOA file for STRING
    goa_file = None
    for path in ["cache/goa_saccharomyces.gaf.gz", "cache/goa_saccharomyces.gaf"]:
        if os.path.exists(path):
            goa_file = path
            print(f"✓ Found GOA file: {path}")
            break
    
    if not goa_file:
        print("⚠ Warning: GOA GAF file not found for STRING mode.")
        print("  To download: wget http://geneontology.org/gene-associations/goa_saccharomyces.gaf.gz -P cache/")
        print("  Or skip STRING mode if GO annotations are not available.")
        string_ppi_exists = False  # Skip STRING if no GO file
    
    print()
    
    success_count = 0
    total_steps = 0
    
    # Step 1: STRING Dataset Pipeline
    if string_ppi_exists:
        total_steps += 1
        cmd = [
            sys.executable, "main.py",
            "--mode", "string",
            "--taxid", "4932",
            "--string-mode", "download",
            "--threshold", "700",
            "--go-file", goa_file,
            "--outdir", "outputs/",
            "--mcl-inflation", "2.0",
            "--alpha", "0.5",
            "--overlap-tau", "0.1",
            "--transfer-tau", "0.0",
            "--lea-population", "30",
            "--lea-evaluations", "500",
            "--lambda-inter", "1.0",
            "--lambda-fragment", "0.5",
            "--random-seed", "42"
        ]
        
        if run_command(cmd, "Step 1: STRING Dataset Pipeline"):
            success_count += 1
    else:
        print("\n⚠ Skipping STRING dataset (files not found)")
    
    # Step 2: Gavin Dataset Pipeline
    if gavin_ppi_exists and go_file_exists:
        total_steps += 1
        cmd = [
            sys.executable, "main.py",
            "--mode", "gavin",
            "--ppi", "gavin2006_socioaffinities_rescaled.txt",
            "--go-file", "GO.txt",
            "--go-use-symbol",
            "--go-taxid", "559292",
            "--outdir", "outputs_gavin/",
            "--mcl-inflation", "2.0",
            "--alpha", "0.5",
            "--overlap-tau", "0.1",
            "--transfer-tau", "0.0",
            "--lea-population", "30",
            "--lea-evaluations", "500",
            "--lambda-inter", "1.0",
            "--lambda-fragment", "0.5",
            "--random-seed", "42"
        ]
        
        if run_command(cmd, "Step 2: Gavin Dataset Pipeline"):
            success_count += 1
    else:
        print("\n⚠ Skipping Gavin dataset (files not found)")
    
    # Step 3: Method Comparison
    total_steps += 1
    cmd = [
        sys.executable, "compare_methods.py",
        "--lea-evaluations", "500"
    ]
    
    if run_command(cmd, "Step 3: Method Comparison"):
        success_count += 1
    
    # Summary
    print("\n" + "="*60)
    print("Summary")
    print("="*60)
    print(f"Completed: {success_count}/{total_steps} steps")
    print()
    
    # List output files
    print("Generated Output Files:")
    print("-"*60)
    
    if os.path.exists("outputs"):
        print("\nSTRING Dataset Outputs (outputs/):")
        for f in sorted(Path("outputs").glob("*.csv")):
            size = f.stat().st_size
            print(f"  {f.name} ({size:,} bytes)")
    
    if os.path.exists("outputs_gavin"):
        print("\nGavin Dataset Outputs (outputs_gavin/):")
        for f in sorted(Path("outputs_gavin").glob("*.csv")):
            size = f.stat().st_size
            print(f"  {f.name} ({size:,} bytes)")
    
    if os.path.exists("community_detection_comparison.csv"):
        size = Path("community_detection_comparison.csv").stat().st_size
        print(f"\nComparison Results:")
        print(f"  community_detection_comparison.csv ({size:,} bytes)")
    
    print()
    print("="*60)
    if success_count == total_steps:
        print("✓ All pipelines completed successfully!")
    else:
        print(f"⚠ Completed {success_count}/{total_steps} steps")
    print("="*60)

if __name__ == "__main__":
    main()

