# Causal ML for Organic Acceptor Molecular Design

**Author:** Samuel Bizimana | JUNIA ISEN  
**Supervisor:** Dr. Kekeli N'KONOU  
**Target Journal:** Nature Machine Intelligence  
**GitHub:** https://github.com/Bsamuel-tech/opv-causal-ml

## What This Project Does

Applies causal machine learning to 15,529 organic acceptor molecules to
identify what structurally causes changes in HOMO energy, LUMO energy,
and optical bandgap. Treatment variable: Hammett-weighted EWG score.
The pipeline goes from dataset construction to quantum chemistry validation
in a fully reproducible workflow.

## Key Results

- EWG score causally lowers LUMO energy by 0.025 eV (DML, p=0.0008)
- IV estimation: F-statistic=304, Wu-Hausman p=0.00175 (IV necessary)
- Note: halogen_count dropped as instrument — it is a component of ewg_weighted
- nonarom_cc is the sole valid instrument (just-identified model)
- Scaffold generalisation gap: mean=-0.0017 eV across 20 bootstrap seeds
- Nuisance diagnostics: R2(ml_l)=0.181, R2(ml_m)=0.713
- Causal forest calibration confirmed (p=0.0004)
- Phase 3 design engine: MAE=0.2155 eV, 6/10 candidates validated
- All 10 candidates synthesizable (SAScore < 4.0)
- Note: ewg_count > 0 filter applied — causal analysis uses 598 molecules

## Analysis Restrictions

All causal models (DML, IV, causal forest) apply df[ewg_count > 0],
restricting analysis to 598 molecules with at least one EWG group.
Nearly all CEP molecules (14,930) have zero EWGs by SMARTS definition
and are excluded from causal estimation. The estimand is the ATE among
EWG-containing organic acceptors, not the full 15,529-molecule population.
See results/tables/supplementary_information.csv Section S0 for details.

## Methods

- Double Machine Learning — DoubleML R package
- Instrumental Variable estimation — ivreg R package (nonarom_cc only)
- Sensitivity analysis — sensemakr R package (benchmark: molecular weight)
- Causal forests and CATE maps — grf R package (4,000 trees)
- Molecular processing and SMARTS reactions — RDKit (Python)
- Quantum chemistry validation — xTB 6.7.1 GFN2

## Repository Contents

**scripts/phase1_data**
- step1_build_dataset.py — builds 599-molecule experimental base
- step2_expand_dataset.py — expands to 15,529 molecules using Harvard CEP
- fix1_add_columns.py — adds measurement_type and halogen_count
- fix2_scaffold_split.py — Murcko scaffold split (70/15/15)
- fix3_acceptor_filter.py — verifies CEP molecules are acceptor-like (99.32% pass)
- fix4_bootstrap_leakage.py — bootstrap leakage gap over 20 seeds
- ewg_utils.py — shared EWG feature computation module

**scripts/phase2_causal**
- dml_analysis.R — Double Machine Learning (DoubleML)
- iv_analysis.R — IV estimation (nonarom_cc instrument only)
- sensitivity_analysis.R — sensemakr sensitivity analysis
- causal_forest.R — causal forest CATE maps (grf)
- econml_comparison.py — EconML vs DoubleML robustness check

**scripts/phase3_design**
- step1_counterfactual_design.py — CATE-driven cyano addition via SMARTS
- step2_xtb_validation.py — xTB GFN2 validation (fixed conformer seed=42)

**results/figures** — 4 publication figures at 300 dpi

**results/tables** — all results tables and supplementary information

**data/processed**
- master_acceptor_dataset.csv — 15,529 molecules (4.38 MB)
  GitHub cannot preview files larger than 1MB in the browser.
  To access: click the file then click View raw to download,
  or run git clone to get all files locally.
- acceptor_base.csv — 599 experimental molecules
- experimental_with_cates.csv — 598 molecules with per-molecule CATEs

## Status

- Phase 0 complete — R 4.6.0 + Python 3.10, fully locked
- Phase 1 complete — 15,529 molecules, bootstrap gap=-0.0017 (20 seeds)
- Phase 2 complete — DML, IV, sensitivity, CATE maps, nuisance diagnostics
- Phase 3 complete — counterfactual design, xTB validation MAE=0.2155 eV
- Phase 4 complete — 4 publication figures, 6-section supplementary

## Data Availability

See DATA_AVAILABILITY.md for full data provenance documentation.
See PATHS_NOTICE.md for notes on running scripts on different machines.

## Reproduce

```r
renv::restore()  # restores all R packages
```

```bash
conda activate causal-mol  # Python environment
conda env create -f environment.yml  # recreate Python environment
```