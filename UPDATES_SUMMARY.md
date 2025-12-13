# Codebase Updates Summary

## Overview
This document summarizes all updates made to ensure consistency with the LEAF-PPI paper definitions, theory, and evaluation methodology.

## 1. Permanence Normalization ✅

**File**: `src/permanence.py`

**Changes**:
- Updated `calculate_permanence()` to properly compute E_max using all clusters (not just external neighbors)
- Added normalization to ensure permanence values lie strictly in range [-1, 1]
- Used linear clipping: `permanence_normalized = max(-1.0, min(1.0, permanence))`
- Added clear comments explaining the normalization choice

**Rationale**: Permanence must be bounded to [-1, 1] for theoretical consistency with literature definitions.

## 2. Functional Dependency (FD) Normalization ✅

**File**: `src/membership_overlap.py`

**Changes**:
- Updated `calculate_functional_dependency()` to normalize FD to [-1, 1] range
- Used `tanh()` normalization to handle unbounded TF-IDF scores
- Added comments explaining that TF-IDF captures functional specificity, not statistical variance
- Ensured FD is comparable in scale to permanence before combination

**Rationale**: FD must be normalized to match permanence scale for meaningful combination in membership score.

## 3. Membership Score Computation ✅

**File**: `src/membership_overlap.py`

**Changes**:
- Updated `calculate_membership()` to use normalized permanence and FD
- Ensured α ∈ [0,1] is respected with clipping
- Added assertions to ensure membership stays in [-1, 1] range
- Updated equation: `score = α * permanence_norm + (1 - α) * fd_norm`

**Rationale**: Membership score combines normalized inputs, ensuring theoretical consistency.

## 4. MCL Community Filtering ✅

**Files**: `src/mcl_clustering.py`, `src/community_comparison.py`, `main.py`, `compare_methods.py`

**Changes**:
- Added `min_cluster_size` parameter (default: 10) to `MCLClustering`
- Filters clusters smaller than `min_cluster_size` before reporting
- Updated all MCL instantiations to include this parameter
- Applied filtering in all fallback methods (cdlib, MCLClustering wrapper)

**Rationale**: Excludes very small clusters (< 10 proteins) to match biological interpretability and published standards (e.g., ~46 for Gavin, not hundreds of trivial clusters).

## 5. External GO-Based Evaluation ✅

**File**: `src/evaluation.py`

**Changes**:
- Added `calculate_go_jaccard_similarity()` function
- Implements Jaccard similarity between GO term sets of detected communities and reference complexes
- Computes mean pairwise Jaccard within clusters if no reference available
- This is EXTERNAL evaluation (separate from FD which is internal/optimization-guiding)

**Rationale**: Provides external functional validation metric for evaluation/reporting, independent of optimization.

## 6. Updated Comparison Script ✅

**Files**: `src/community_comparison.py`, `compare_methods.py`

**Changes**:
- Updated `evaluate_communities()` to include `mean_go_jaccard` metric
- Updated `compare_all_methods()` to accept `protein_go_terms` parameter
- Updated all method calls to pass GO terms for external evaluation
- Added `mean_go_jaccard` to all result dictionaries

**Rationale**: Enables comprehensive evaluation including external GO-based metrics.

## 7. CSV Output Format ✅

**Files**: `src/community_comparison.py`, `generate_updated_results.py`

**Changes**:
- CSV columns include: `dataset`, `method`, `overlapping`, `num_nodes`, `num_edges`, `num_communities`, `avg_community_size`, `modularity`, `nmi`, `overlapping_nmi`, `conductance`, `mean_go_jaccard`, `runtime_sec`, `parameters`
- Created `generate_updated_results.py` script to regenerate all results
- Output files: `results_string_updated.csv`, `results_gavin_updated.csv`, `community_detection_comparison.csv`

**Rationale**: Standardized output format with all required metrics.

## 8. Reproducibility ✅

**Files**: All relevant files

**Changes**:
- Fixed random seed: 42 (default)
- Ensured identical preprocessing for all methods
- Added clear inline comments explaining evaluation choices
- All normalization choices documented

**Rationale**: Ensures reproducible results for scientific publication.

## Key Files Modified

1. `src/permanence.py` - Permanence normalization
2. `src/membership_overlap.py` - FD normalization, membership computation
3. `src/mcl_clustering.py` - MCL filtering
4. `src/evaluation.py` - External GO evaluation
5. `src/community_comparison.py` - Updated comparison framework
6. `compare_methods.py` - Updated to pass GO terms
7. `main.py` - Added min_cluster_size parameter
8. `generate_updated_results.py` - New script for regenerating results

## Next Steps

1. Run `python generate_updated_results.py` to regenerate all CSV files
2. Verify all metrics are within expected ranges
3. Check that MCL community counts match expected values (~46 for Gavin)
4. Validate that permanence and FD values are in [-1, 1] range
5. Review generated CSV files for consistency

## Notes

- All changes maintain backward compatibility where possible
- Normalization choices are conservative (linear clipping for permanence, tanh for FD)
- MCL filtering is configurable via `min_cluster_size` parameter
- External GO evaluation is optional (returns None if GO terms unavailable)

