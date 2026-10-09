# Pipeline Reference

Complete record of all scripts, their inputs, outputs, and row counts as of 2026-10-09.

## Phase 1: Data Preparation

### fix1_feature_engineering.py
- **Input:** raw merged CSV (experimental literature + CEP download)
- **Output:** `data/interim/fix1_features.csv`
- **Creates columns:** `measurement_type` ("experiment" or "DFT_B3LYP"), `halogen_count`, `nonarom_cc`, `ewg_count`, `ewg_weighted`, `mol_weight`, `n_arom_rings`, `conj_length`, `logp`
- **Note:** `measurement_type` value is "experiment" not "experimental"

### fix2_scaffold_split.py
- **Input:** `data/interim/fix1_features.csv`
- **Output:** `data/interim/fix2_with_scaffolds.csv`
- **Creates columns:** `scaffold` (Murcko), `split` (train/val/test)
- **Guarantee:** no scaffold appears in more than one split

### fix3_acceptor_filter.py
- **Input:** `data/interim/fix2_with_scaffolds.csv`
- **Output:** `data/interim/fix3_acceptors_only.csv`
- **Rows:** 15,403 (CEP molecules filtered to LUMO < -2.0 eV, bandgap < 3.5 eV, HOMO < -4.5 eV)
- **Filter report:** `results/tables/cep_acceptor_filter_report.csv`

### fix4_bootstrap_leakage.py
- **Input:** `data/interim/fix3_acceptors_only.csv`
- **Outputs:**
  - `data/processed/master_acceptor_dataset.csv` (15,403 rows — canonical master for all downstream)
  - `results/tables/bootstrap_leakage_gap.csv` (per-seed gaps, both populations)
  - `results/tables/bootstrap_leakage_summary.csv` (labelled summary, two sections)
- **Section A:** full dataset (~15k) bootstrap — informative for prediction model only
- **Section B:** causal subset (598 molecules, measurement_type=="experiment" & ewg_count>0) — relevant for DML validity

**Note on iv_evaluation_table4.csv:** This file was removed (commit 292f723). Git history shows it first appeared in commit b74982b with statistics referencing two instruments including halogen_count. halogen_count is a component of ewg_weighted and was never a valid instrument. The numbers in that file cannot be traced to a reproducible script run and were excluded from all reported results.

## Phase 2: Causal Inference

All Phase 2 scripts read from `data/processed/master_acceptor_dataset.csv` and filter to the causal estimation subset: `measurement_type == "experiment" & ewg_count > 0` → **598 rows, 288 unique scaffolds**.

### dml_analysis.R
- **Fold assignment:** cyclic round-robin after random permutation (replace=FALSE). Replaces earlier sample(1:5, replace=TRUE) which risked empty folds.
- **Guards:** aborts on zero rows, NA scaffolds, scaffold spanning two folds, empty fold, train/test overlap
- **Output:** `results/tables/dml_ewg_corrected.csv`
  - Columns: treatment, outcome, ATE, SE, p_value, CI_lower, CI_upper, fold_scheme, n_obs, n_scaffolds, n_folds
- **Diagnostic output:** `results/tables/dml_fold_diagnostics.csv` (per-fold molecule and scaffold counts)

### iv_analysis.R
- **Instrument:** `nonarom_cc` only (just-identified; Sargan not applicable)
- **halogen_count** is a confounder, not an instrument — it is algebraically a component of ewg_weighted
- **Output:** `results/tables/iv_results.csv`

### sensitivity_analysis.R
- **Method:** sensemakr (Cinelli & Hazlett 2020)
- **Benchmark covariate:** mol_weight
- **Output:** `results/tables/sensemakr_results.csv`, `results/figures/sensitivity_contour.png`

### causal_forest.R
- **Method:** grf package, 4000 trees
- **Output:** `results/tables/cate_results.csv`, `results/figures/cate_heatmap.png`

## Phase 3: Counterfactual Design and Validation

### counterfactual_generator.py
- **Input:** causal forest CATE predictions + top candidate SMILES
- **Method:** SMARTS reaction to add cyano groups; RDKit SAScore filter (< 4.0)
- **Output:** `results/tables/counterfactual_candidates.csv`

### xtb_validation.py
- **Input:** `counterfactual_candidates.csv`
- **Method:** ETKDGv3 conformer generation + xTB GFN2 single-point
- **Output:** `results/tables/xtb_validation_table.csv`
- **Result:** MAE = 0.219 eV; 5/10 candidates pass SAScore < 4.0

## Key Numbers (as of audit 2026-10-09)

| Quantity | Value |
|----------|-------|
| Total molecules (master) | 15,403 |
| Causal estimation subset | 598 molecules |
| Unique scaffolds (causal) | 288 |
| DML HOMO ATE | -0.0179 eV, p=0.059 |
| DML LUMO ATE | -0.0161 eV, p=0.051 |
| DML Bandgap ATE | -0.0025 eV, p=0.798 |
| xTB validation MAE | 0.219 eV |
| Candidates passing SAScore | 5/10 |
| Bootstrap gap mean (full) | 0.0122 |