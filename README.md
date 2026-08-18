# Causal ML for Organic Acceptor Molecular Design

**Author:** Samuel Bizimana | JUNIA ISEN  
**Supervisor:** Dr. Kekeli N'KONOU  
**Target Journal:** Nature Machine Intelligence  

## What This Project Does

Applies causal machine learning to 15,529 organic acceptor molecules to
identify what structurally causes changes in HOMO energy. Treatment variable is the Hammett-weighted EWG score.
The pipeline runs from dataset construction to quantum chemistry validation
in a fully reproducible workflow.

## Key Results

EWG score causally lowers LUMO energy by 0.025 eV (DML, p=0.0008).
IV estimation with nonarom_cc as sole valid instrument gives F-statistic=216.
halogen_count was dropped as instrument because it is algebraically a component
of ewg_weighted and violates the exclusion restriction. Wu-Hausman p=0.433
is inconclusive at this sample size and is reported honestly as a limitation.
Bootstrap leakage gap mean=-0.0017 across 20 seeds (95% CI [-0.046, +0.044]),
all below the 0.15 threshold, confirming no data leakage.
Nuisance diagnostics: R2(ml_l)=0.181 RMSE=0.133 eV, R2(ml_m)=0.713.
Causal forest calibration confirmed at p=0.0004.
Phase 3 design engine achieves MAE=0.2155 eV with 6/10 candidates validated
and all 10 synthesizable with SAScore below 4.0.

## Analysis Restriction

All causal models apply ewg_count > 0, restricting causal estimation to
598 molecules with at least one EWG group. Nearly all 14,930 CEP molecules
have zero EWGs by SMARTS definition and contribute to scaffold diversity
but not to causal estimation. The estimand is the ATE among EWG-containing
organic acceptors, not the full 15,529-molecule population.
See results/tables/supplementary_information.csv Section S0 for full details.

## Methods

Double Machine Learning via DoubleML in R.
Instrumental Variable estimation via ivreg in R using nonarom_cc only.
Sensitivity analysis via sensemakr in R with molecular weight as benchmark.
Causal forests and CATE maps via grf in R with 4,000 trees.
Molecular processing and SMARTS reactions via RDKit in Python.
Quantum chemistry validation via xTB 6.7.1 GFN2 with fixed conformer seed=42.

## Repository Contents

scripts/phase1_data contains step1_build_dataset.py, step2_expand_dataset.py,
fix1_add_columns.py, fix2_scaffold_split.py, fix3_acceptor_filter.py,
fix4_bootstrap_leakage.py, and ewg_utils.py as a shared EWG module.

scripts/phase2_causal contains dml_analysis.R, iv_analysis.R,
sensitivity_analysis.R, causal_forest.R, and econml_comparison.py.

scripts/phase3_design contains step1_counterfactual_design.py and
step2_xtb_validation.py with fixed conformer seed for reproducibility.

results/figures contains four publication figures at 300 dpi.
results/tables contains all results tables and supplementary information.

data/processed contains master_acceptor_dataset.csv with 15,529 molecules
at 4.38 MB. GitHub cannot preview files larger than 1MB in the browser.
To access the file click it on GitHub then click View raw to download,
or run git clone to get all files locally.

## Status

Phase 0 complete: R 4.6.0 + Python 3.10, all packages locked in
renv.lock, requirements.txt, and environment.yml.
Phase 1 complete: 15,529 molecules, bootstrap leakage gap mean=-0.0017
across 20 seeds, CEP acceptor filter 99.32% pass rate.
Phase 2 complete: DML, IV (corrected instrument), sensitivity analysis,
CATE maps, nuisance diagnostics, EconML robustness check.
Phase 3 complete: counterfactual design engine using per-molecule CATE,
xTB validation MAE=0.2155 eV, 6/10 candidates validated, 10/10 synthesizable.
Phase 4 complete: four publication figures, seven-section supplementary.

## Data Availability

See DATA_AVAILABILITY.md for full data provenance documentation including
CEP database URL, exact filters applied, and experimental dataset description.
See PATHS_NOTICE.md for notes on running scripts on different machines.

## Reproduce

```r
renv::restore()
```

```bash
conda activate causal-mol
conda env create -f environment.yml
```
## What This Project Does

Applies causal machine learning to 15,529 organic acceptor molecules to
identify what structurally causes changes in HOMO energy, LUMO energy,
and optical bandgap. Treatment variable is the Hammett-weighted EWG score.
The pipeline runs from dataset construction to quantum chemistry validation
in a fully reproducible workflow.

## Key Results

EWG score causally lowers LUMO energy by 0.025 eV (DML, p=0.0008).
IV estimation gives F-statistic=304 and Wu-Hausman p=0.00175, confirming
IV was necessary. The sole valid instrument is nonarom_cc since halogen_count
is a component of the treatment variable ewg_weighted and was dropped.
Bootstrap leakage gap is mean=0.0017 across 20 seeds, confirming no leakage.
Nuisance diagnostics: R2(ml_l)=0.181, R2(ml_m)=0.713.
Causal forest calibration confirmed at p=0.0004.
Phase 3 design engine achieves MAE=0.2155 eV with 6/10 candidates validated
and all 10 synthesizable with SAScore below 4.0.

## Analysis Restriction

All causal models apply ewg_count > 0, restricting analysis to 598 molecules
with at least one EWG group. Nearly all CEP molecules have zero EWGs by SMARTS
definition and are excluded from causal estimation. The estimand is the ATE
among EWG-containing organic acceptors, not the full 15,529-molecule population.
See results/tables/supplementary_information.csv Section S0 for full details.

## Methods

Double Machine Learning via DoubleML in R.
Instrumental Variable estimation via ivreg in R using nonarom_cc only.
Sensitivity analysis via sensemakr in R with molecular weight as benchmark.
Causal forests and CATE maps via grf in R with 4,000 trees.
Molecular processing and SMARTS reactions via RDKit in Python.
Quantum chemistry validation via xTB 6.7.1 GFN2 with fixed conformer seed.

## Repository Contents

scripts/phase1_data contains step1_build_dataset.py, step2_expand_dataset.py,
fix1_add_columns.py, fix2_scaffold_split.py, fix3_acceptor_filter.py,
fix4_bootstrap_leakage.py, and ewg_utils.py as a shared EWG module.

scripts/phase2_causal contains dml_analysis.R, iv_analysis.R,
sensitivity_analysis.R, causal_forest.R, and econml_comparison.py.

scripts/phase3_design contains step1_counterfactual_design.py and
step2_xtb_validation.py with fixed conformer seed for reproducibility.

results/figures contains four publication figures at 300 dpi.
results/tables contains all results tables and supplementary information.

data/processed contains master_acceptor_dataset.csv with 15,529 molecules
at 4.38 MB. GitHub cannot preview files larger than 1MB in the browser.
To access the file, click it on GitHub then click View raw to download,
or run git clone to get all files locally. Also contains acceptor_base.csv
with 599 experimental molecules and experimental_with_cates.csv with
598 molecules including per-molecule CATEs from the causal forest.

## Status

Phase 0 complete: R 4.6.0 + Python 3.10, all packages locked.
Phase 1 complete: 15,529 molecules, bootstrap leakage gap mean=0.0017.
Phase 2 complete: DML, IV, sensitivity, CATE maps, nuisance diagnostics.
Phase 3 complete: counterfactual design engine, xTB validation MAE=0.2155 eV.
Phase 4 complete: four publication figures, six-section supplementary.

## Data Availability

See DATA_AVAILABILITY.md for full data provenance documentation.
See PATHS_NOTICE.md for notes on running scripts on different machines.

## Reproduce

```r
renv::restore()
```

```bash
conda activate causal-mol
conda env create -f environment.yml
```
