# PIPELINE.md: OPV Causal ML Pipeline
Run scripts in the order listed. Each step reads from the previous step's output.

## Phase 1: Data

| Step | Script | Input | Output | Rows (approx) |
|------|--------|-------|--------|---------------|
| 1 | `scripts/phase1_data/step1_build_dataset.py` | `data/raw/Data_Merged_with_SMILES.xlsx` | `data/interim/step1_experimental_raw.csv` | 1,573 |
| 2 | `scripts/phase1_data/step2_expand_dataset.py` | `data/interim/step1_experimental_raw.csv` + `data/raw/cep_sample_15000.csv` | `data/interim/step2_merged_raw.csv` | ~14,000 |
| 3 | `scripts/phase1_data/fix1_add_columns.py` | `data/interim/step2_merged_raw.csv` | `data/interim/fix1_with_columns.csv` | ~14,000 |
| 4 | `scripts/phase1_data/fix2_scaffold_split.py` | `data/interim/fix1_with_columns.csv` | `data/interim/fix2_with_scaffolds.csv` | ~14,000 |
| 5 | `scripts/phase1_data/fix3_acceptor_filter.py` | `data/interim/fix2_with_scaffolds.csv` | `data/interim/fix3_acceptors_only.csv` | 15,403 |
| 6 | `scripts/phase1_data/fix4_bootstrap_leakage.py` | `data/interim/fix3_acceptors_only.csv` | `data/processed/master_acceptor_dataset.csv` + `results/tables/bootstrap_leakage_gap.csv` + `results/tables/bootstrap_leakage_summary.csv` | 15,403 |

## Phase 2: Causal Inference

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 7 | `scripts/phase2_causal/dml_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/dml_ewg_corrected.csv`, `results/tables/dml_fold_diagnostics.csv` |
| 8 | `scripts/phase2_causal/causal_forest.R` | `data/processed/master_acceptor_dataset.csv` | `data/processed/experimental_with_cates.csv`, `results/tables/cate_results_corrected.csv`, `results/figures/cate_ewg_homo_distribution.png`, `results/figures/cate_ewg_vs_molweight.png` |
| 9 | `scripts/phase2_causal/iv_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/iv_results.csv` |
| 10 | `scripts/phase2_causal/sensitivity_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/sensitivity_results.csv`, `results/figures/sensitivity_contour.png` |

## Phase 3: Counterfactual Design

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 11 | `scripts/phase3_design/step1_counterfactual_design.py` | `data/processed/experimental_with_cates.csv` | `results/tables/phase3_top10_candidates.csv` |
| 12 | `scripts/phase3_design/step2_xtb_validation.py` | `results/tables/phase3_top10_candidates.csv` | `results/tables/phase3_xtb_validation_final.csv` |

## Key results (verified from R run 2026-10-09)

| Outcome | ATE (eV) | SE | p-value | CI | Significant |
|---------|----------|----|---------|----|-------------|
| HOMO | -0.02223 | 0.00861 | 0.0098 | [-0.0391, -0.0054] | Yes |
| LUMO | -0.01538 | 0.00867 | 0.0759 | [-0.0324, +0.0016] | No |
| Bandgap | +0.00173 | 0.00974 | 0.8593 | [-0.0174, +0.0208] | No |

n_obs=598 | n_scaffolds=288 | fold_scheme=scaffold-grouped | n_folds=5

## Raw data notes

- `data/raw/Data_Merged_with_SMILES.xlsx`: 1,573 experimental acceptor molecules
- `data/raw/cep_sample_15000.csv`: 15,000-row reproducible sample from Harvard CEPDB, sampled with `random_state=42` / `set.seed(42)`

## Known issues and decisions

- `measurement_type` column values are `"experiment"` and `"DFT_B3LYP"` (not `"experimental"`). All Phase 2 R scripts filter on `"experiment"`.
- `iv_evaluation_table4.csv` has been removed (commit 292f723). The IV results table is `results/tables/iv_results.csv`, produced by `iv_analysis.R` with `nonarom_cc` as the sole instrument. `halogen_count` is not used as an instrument (exclusion restriction violation — it is a component of `ewg_weighted`). The just-identified model has no Sargan test.
- `causal_forest.R` produces two figures: `cate_ewg_homo_distribution.png` and `cate_ewg_vs_molweight.png`. There is no `cate_heatmap.png`.
- DML cross-fitting uses scaffold-grouped folds (cyclic round-robin assignment, not random) to prevent train/test scaffold overlap. Fold sizes: 183 / 91 / 106 / 102 / 116 molecules.
- Causal subset (598 molecules) is distinct from the full dataset (15,403 molecules). All Phase 2 causal estimates apply only to the causal subset.
- xTB must be installed separately and available on PATH for `step2_xtb_validation.py` to run.