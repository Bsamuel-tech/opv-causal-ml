# Causal Machine Learning for Organic Acceptor Design

**Author:** Samuel Bizimana | JUNIA ISEN  
**Supervisor:** Dr. Kekeli N'KONOU  
**Target Journal:** Nature Machine Intelligence  
**GitHub:** https://github.com/Bsamuel-tech/opv-causal-ml

**Pipeline:** https://github.com/Bsamuel-tech/opv-causal-ml/blob/main/PIPELINE.md

## What This Project Does

Applies causal machine learning to 15,403 organic acceptor molecules to
identify what structurally causes changes in HOMO energy, LUMO energy,
and optical bandgap. Treatment variable is the Hammett-weighted EWG score.
The pipeline runs from dataset construction to quantum chemistry validation
in a fully reproducible workflow.

## Key Results

EWG score causally lowers HOMO energy by -0.0222 eV (DML, p=0.010;
95% CI [-0.039, -0.005] excludes zero — statistically significant).
LUMO shows a negative but inconclusive estimate (-0.0154 eV, p=0.076;
CI [-0.032, +0.002] crosses zero). Bandgap shows no effect (+0.0017 eV,
p=0.859). All three estimates come from scaffold-grouped 5-fold
cross-fitting on 598 experimental molecules across 288 unique scaffolds,
using cyclic round-robin fold assignment (verified run 2026-10-09).
No valid instrument was identified for IV estimation. halogen_count is
retained as a confounder because it contains structural information beyond
its weighted contribution to the Hammett EWG score. nonarom_cc was tested
as an instrument but its exclusion restriction is chemically uncertain —
non-aromatic C=C bonds may directly affect conjugation length, which
determines orbital energies. IV results are inconclusive; DML is the
primary estimate.
Bootstrap leakage gap mean=0.0122 (full dataset, 20 seeds, all below 0.15
threshold). The causal estimation subset (598 molecules) is bootstrapped
separately — see Section B of bootstrap_leakage_summary.csv.
Nuisance diagnostics: R2(ml_l)=0.181 RMSE=0.133 eV, R2(ml_m)=0.713.
Causal forest calibration confirmed (p=0.0004). Heterogeneity not confirmed
(p=0.295) — reported as an honest negative finding.
Phase 3 design engine achieves MAE=0.219 eV; 5/10 candidates pass both
xTB validation and SAScore < 4.0.

## Analysis Restriction

All causal models apply ewg_count > 0, restricting causal estimation to
598 molecules with at least one EWG group. Nearly all CEP molecules
have zero EWGs by SMARTS definition and contribute to scaffold diversity
but not to causal estimation. The estimand is the ATE among EWG-containing
organic acceptors, not the full 15,403-molecule population.
See results/tables/supplementary_information.csv Section S0 for full details.

## Methods

Double Machine Learning via DoubleML in R. Scaffold-grouped 5-fold
cross-fitting with cyclic round-robin fold assignment (guarantees no empty
folds). Filter: measurement_type == "experiment" AND ewg_count > 0.
IV estimation attempted but no valid instrument identified — reported as
a limitation. Sensitivity analysis via sensemakr in R with molecular weight
as benchmark. Causal forests and CATE maps via grf in R with 4,000 trees.
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

results/figures contains four publication figures at 300 dpi, including
cate_ewg_homo_distribution.png and cate_ewg_vs_molweight.png.
results/tables contains all results tables and supplementary information,
including dml_ewg_corrected.csv and dml_fold_diagnostics.csv.
results/models contains causal_forest_ewg_homo.rds (5MB trained model).

data/processed contains master_acceptor_dataset.csv with 15,403 molecules.
GitHub cannot preview files larger than 1MB in the browser.
To access the file click it on GitHub then click View raw to download,
or run git clone to get all files locally.

## Status

Phase 0 complete: R 4.6.0 + Python 3.10, all packages locked in
renv.lock, requirements.txt, and environment.yml.
Phase 1 complete: 15,403 molecules, bootstrap leakage gap mean=0.0122
(full dataset, 20 seeds), CEP acceptor filter 99.32% pass rate.
Phase 2 complete: DML verified run 2026-10-09 (scaffold-grouped cyclic
folds, 598 molecules, 288 scaffolds). HOMO effect significant (p=0.010),
LUMO inconclusive (p=0.076), Bandgap no effect (p=0.859). Fold diagnostics
in dml_fold_diagnostics.csv. IV no valid instrument — reported as
limitation. Sensitivity analysis and CATE maps complete.
Phase 3 complete: counterfactual design using per-molecule CATE,
xTB validation MAE=0.219 eV, 5/10 pass SAScore < 4.0.
Phase 4 complete: four publication figures, seven-section supplementary.

## Data Availability

See DATA_AVAILABILITY.md for full data provenance documentation.
See PATHS_NOTICE.md for notes on running scripts on different machines.

## Reproduce

```r
renv::restore()
source("scripts/phase2_causal/dml_analysis.R")
```

```bash
conda activate causal-mol
conda env create -f environment.yml
python scripts/phase1_data/fix4_bootstrap_leakage.py
```