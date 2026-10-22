# PIPELINE.md: OPV Causal ML Pipeline
Run scripts in the order listed. Each step reads from the previous step's output.

## Phase 1: Data

| Step | Script | Input | Output | Rows |
|------|--------|-------|--------|------|
| 1 | `scripts/phase1_data/step1_build_dataset.py` | `data/raw/Data_Merged_with_SMILES.xlsx` | `data/interim/step1_experimental_raw.csv` | 599 (output rows; input xlsx has 1,573 including header/duplicates) |
| 2 | `scripts/phase1_data/step2_expand_dataset.py` | `data/interim/step1_experimental_raw.csv` + `data/raw/cep_sample_15000.csv` | `data/interim/step2_merged_raw.csv` | 15,528 |
| 3 | `scripts/phase1_data/fix1_add_columns.py` | `data/interim/step2_merged_raw.csv` | `data/interim/fix1_with_columns.csv` | 15,528 |
| 4 | `scripts/phase1_data/fix2_scaffold_split.py` | `data/interim/fix1_with_columns.csv` | `data/interim/fix2_with_scaffolds.csv` | 15,528 |
| 5 | `scripts/phase1_data/fix3_acceptor_filter.py` | `data/interim/fix2_with_scaffolds.csv` | `data/interim/fix3_acceptors_only.csv` | 15,403 |
| 6 | `scripts/phase1_data/fix4_bootstrap_leakage.py` | `data/interim/fix3_acceptors_only.csv` | `data/processed/master_acceptor_dataset.csv` + `results/tables/bootstrap_leakage_gap.csv` + `results/tables/bootstrap_leakage_summary.csv` | 15,403 |

## Phase 2: Causal Inference

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 7 | `scripts/phase2_causal/dml_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/dml_ewg_corrected.csv`, `results/tables/dml_fold_diagnostics.csv`, `results/tables/dml_cluster_robust.csv` |
| 8 | `scripts/phase2_causal/causal_forest.R` | `data/processed/master_acceptor_dataset.csv` | `data/processed/experimental_with_cates.csv`, `results/tables/cate_results_corrected.csv`, `results/figures/cate_ewg_homo_distribution.png`, `results/figures/cate_ewg_vs_molweight.png` |
| 9 | `scripts/phase2_causal/iv_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/iv_results.csv`, `results/tables/iv_diagnostics.csv` |
| 10 | `scripts/phase2_causal/sensitivity_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/sensitivity_results.csv`, `results/figures/sensitivity_contour.png` |

## Phase 3: Counterfactual Design

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 11 | `scripts/phase3_design/step1_counterfactual_design.py` | `data/processed/experimental_with_cates.csv` | `results/tables/phase3_top10_candidates.csv` |
| 12 | `scripts/phase3_design/step2_xtb_validation.py` | `results/tables/phase3_top10_candidates.csv` | `results/tables/phase3_xtb_validation_final.csv` |

## Key results (from committed code, run 2026-10-19, n_rep=20 scaffold-grouped; not independently reproduced)

| Outcome | ATE (eV) | SE IID / SE CR | p IID / p CR | CI (IID) | Inference |
|---------|----------|----------------|--------------|----------|-----------|
| HOMO | -0.01828 | 0.00899 / 0.00940 | 0.0420 / 0.0520 | [-0.03589, -0.00066] | Borderline |
| LUMO | -0.01633 | 0.00830 / 0.00837 | 0.0491 / 0.0510 | [-0.03260, -0.00006] | Borderline |
| Bandgap | -0.00228 | 0.00969 / 0.01139 | 0.8138 / 0.8410 | [-0.02128, +0.01671] | No effect |

IID = DoubleML i.i.d. SE; CR = post-hoc scaffold cluster-robust sandwich SE (Cameron & Miller 2015, rowMeans across 20 reps).
Cluster-robust CIs: HOMO [-0.03671, +0.00016]; LUMO [-0.03273, +0.00007] — both marginally cross zero.
Bonferroni-adjusted p-values (3 outcomes): HOMO ~0.13; LUMO ~0.15; both non-significant after correction.

n_obs=598 | n_scaffolds=288 | fold_scheme=scaffold-grouped-cyclic | n_folds=5 | n_rep=20

Causal forest calibration: mean effect confirmed (p=0.0039). Heterogeneity not confirmed (differential.forest.prediction p=0.629, one-sided HC3, from test_calibration(cf)).

## Lineage manifest

| File | MD5 | Rows | Cols | Description |
|------|-----|------|------|-------------|
| `data/raw/Data_Merged_with_SMILES.xlsx` | `99b06236592b7663fd3272cb2608d640` | 1,573 | — | Raw experimental acceptors including duplicates |
| `data/interim/step1_experimental_raw.csv` | `cb04b2a9a2665daa557f36fec26f8ead` | 599 | 12 | Deduplicated experimental acceptors |
| `data/interim/step2_merged_raw.csv` | — | 15,528 | 12 | Merged experimental + CEP before column enrichment |
| `data/interim/fix1_with_columns.csv` | — | 15,528 | 15 | After column standardization |
| `data/interim/fix2_with_scaffolds.csv` | — | 15,528 | 17 | After scaffold assignment |
| `data/interim/fix3_acceptors_only.csv` | `664009c84abab43ee3c3a6fc0600c419` | 15,403 | 17 | Full merged dataset after acceptor filter |
| `data/processed/master_acceptor_dataset.csv` | `664009c84abab43ee3c3a6fc0600c419` | 15,403 | 17 | Final master dataset (copy of fix3 after leakage check) |
| `results/tables/dml_ewg_corrected.csv` | `35979fdc752c93f85fc7efa4a73a26e5` | 3 | 15 | DML ATE results (one row per outcome), n_rep=20 |
| `results/tables/dml_cluster_robust.csv` | `f7b21e1dee183851579d7c5255e57420` | 3 | 8 | Scaffold cluster-robust SEs (outcome, ATE, SE_iid, SE_cluster, p_cluster, CI_lower_cr, CI_upper_cr, n_clusters) |
| `results/tables/dml_fold_diagnostics.csv` | `5af92d85c6a91bb3273729f56f4f8caa` | 5 | 3 | Fold size diagnostics for seed 1 |
| `results/tables/iv_results.csv` | `99464bf1d97b44ff76f7a3897222601e` | 2 | 4 | IV coefficients (one row per outcome) |
| `results/tables/iv_diagnostics.csv` | `cf7d54527988bd7166950fdfe3a0922b` | 2 | 5 | IV diagnostics: F-stat, Wu-Hausman p (one row per outcome) |
| `results/tables/cate_results_corrected.csv` | `f45d5dc5cd6c18801265443b6d8a9387` | 598 | 7 | Per-molecule CATE estimates with scaffold column |
| `data/processed/experimental_with_cates.csv` | `f45d5dc5cd6c18801265443b6d8a9387` | 598 | 7 | Causal subset with CATE column attached |
| `results/tables/bootstrap_leakage_summary.csv` | `71c5ff401ecb5d294a6427d6b0e665f7` | 17 | 2 | Section A and B leakage gap summary |
| `results/tables/bootstrap_leakage_gap.csv` | `d1fcaa509bf73d7ca98f900fa962643e` | 20 | 7 | Per-seed leakage gaps (full and causal subset) |
| `results/tables/phase3_xtb_validation_final.csv` | `2a6ef1b58b242e40638af276432e3bdc` | 10 | 14 | xTB-validated counterfactual candidates |
| `results/tables/supplementary_information.csv` | `ee4e82b26e39104602712bfff5ae687b` | — | — | All supplementary sections S0-S6 |

Note: MD5s marked — not yet verified for interim files (steps 2-4). Run `scripts/generate_manifest_md5.sh` after any pipeline rerun to refresh hashes.

## Raw data notes

- `data/raw/Data_Merged_with_SMILES.xlsx`: 1,573 experimental acceptor molecules
- `data/raw/cep_sample_15000.csv`: 15,000-row reproducible sample from Harvard CEPDB, sampled with `random_state=42` / `set.seed(42)`

## Known issues and decisions

- `measurement_type` column values are `"experiment"` and `"DFT_B3LYP"` (not `"experimental"`). All Phase 2 R scripts filter on `"experiment"`.
- `iv_evaluation_table4.csv` has been removed (commit 292f723). The IV results table is `results/tables/iv_results.csv`, produced by `iv_analysis.R` with `nonarom_cc` as the sole instrument. `halogen_count` is not used as an instrument (exclusion restriction violation — it is a component of `ewg_weighted`). The just-identified model has no Sargan test.
- `causal_forest.R` produces two figures: `cate_ewg_homo_distribution.png` and `cate_ewg_vs_molweight.png`. There is no `cate_heatmap.png`.
- DML cross-fitting uses scaffold-grouped folds (cyclic round-robin assignment, not random) to prevent train/test scaffold overlap. Fold sizes at seed 1: 118 / 131 / 151 / 79 / 119 molecules. Sizes vary across seeds by design.
- `causal_forest.R` passes `clusters = scaffold_cluster_ids` (288 integer cluster IDs) to `grf::causal_forest()` for cluster-robust variance estimation.
- `dml_analysis.R` applies a post-hoc scaffold cluster-robust sandwich SE correction (Cameron & Miller 2015, rowMeans across 20 reps). Output in `results/tables/dml_cluster_robust.csv`. The i.i.d. SEs in `dml_ewg_corrected.csv` are also retained for comparison.
- Section B bootstrap (causal subset, n=598): gap mean=0.316, 95% CI [0.164, 0.486]. Large gap confirms scaffold leakage in naive random splits, validating scaffold-grouped cross-fitting in all Phase 2 models.
- Causal subset (598 molecules) is distinct from the full dataset (15,403 molecules). All Phase 2 causal estimates apply only to the causal subset.
- xTB must be installed separately and available on PATH for `step2_xtb_validation.py` to run.
- IV exclusion restriction for `nonarom_cc` is chemically uncertain. Non-aromatic C=C bonds may directly affect conjugation length and orbital energies. IV coefficients (~-0.066 eV) are approximately 3.7x larger than DML estimates; IV supports the direction but not the magnitude. IV results reported for transparency; DML is the primary estimate.
- Wu-Hausman p<0.05 for HOMO and LUMO confirms endogeneity (OLS and IV inconsistent). This supports use of IV estimation but does not validate the instrument itself.