# Causal ML for Organic Acceptor Molecular Design
## Overview

Causal machine learning pipeline for organic acceptor molecules.
Treatment: Hammett-weighted EWG score.
Outcomes: HOMO energy, LUMO energy, optical bandgap.
Dataset: 15,529 molecules. Validation: xTB GFN2 quantum chemistry.

## Repository Contents

**scripts/phase1_data**
- step1_build_dataset.py — builds the 599-molecule experimental base
- step2_expand_dataset.py — expands to 15,529 molecules using Harvard CEP
- fix1_add_columns.py — adds measurement type and halogen count columns
- fix2_scaffold_split.py — Murcko scaffold split (70/15/15)

**scripts/phase2_causal**
- dml_analysis.R — Double Machine Learning (DoubleML, R)
- iv_analysis.R — Instrumental Variable estimation (ivreg, R)
- econml_comparison.py — EconML vs DoubleML robustness check (Python)

**scripts/phase3_design**
- step1_counterfactual_design.py — CATE-driven cyano addition via SMARTS
- step2_xtb_validation.py — xTB GFN2 validation of top 10 candidates

**results/figures**
- Figure1_dataset.png — dataset composition and scaffold split
- Figure2_cate_maps.png — CATE heatmap and scatter with confidence bands
- fig_phase3_design_loop.png — counterfactual design loop diagram
- fig_phase3_scatter.png — ML predicted vs xTB computed HOMO

**results/tables**
- dml_ewg_corrected.csv — DML results for HOMO, LUMO, bandgap
- iv_diagnostics_complete.csv — F-stat, Sargan, Wu-Hausman results
- econml_vs_doubleml_comparison.csv — cross-framework robustness check
- phase3_xtb_validation_final.csv — xTB validation of 10 candidates
- supplementary_information.csv — 6-section supplementary document

**data/processed**
- acceptor_base.csv — 599 experimental acceptor molecules
- master_acceptor_dataset.csv — full 15,529 molecule dataset

## Key Results

- EWG causally lowers LUMO by 0.025 eV (p = 0.0008)
- IV F-statistic = 304, Wu-Hausman p = 0.00175 (IV necessary)
- Sensemakr robustness value = 0.217
- Scaffold generalisation gap = 0.050 (no data leakage)
- Phase 3 MAE = 0.220 eV, 5/10 candidates validated, 10/10 synthesizable

## Reproduce

```r
renv::restore()
```

```bash
conda activate causal-mol
```

## Note

Large files (raw dataset 563 MB, CEP data 563 MB) are kept locally only.
