# PIPELINE.md: OPV Causal ML Pipeline
Run scripts in the order listed. Each step reads from the previous step's output.

## Phase 1: Data

| Step | Script | Input | Output | Rows (approx) |
|------|--------|-------|--------|---------------|
| 1 | `scripts/phase1_data/step1_build_dataset.py` | `data/raw/Data_Merged_with_SMILES.xlsx` | `data/interim/step1_experimental_raw.csv` | 1,573 |
| 2 | `scripts/phase1_data/step2_expand_dataset.py` | `data/interim/step1_experimental_raw.csv` + `data/raw/cep_sample_15000.csv` | `data/interim/step2_merged_raw.csv` | ~14,000 |
| 3 | `scripts/phase1_data/fix1_add_columns.py` | `data/interim/step2_merged_raw.csv` | `data/interim/fix1_with_columns.csv` | ~14,000 |
| 4 | `scripts/phase1_data/fix2_scaffold_split.py` | `data/interim/fix1_with_columns.csv` | `data/interim/fix2_with_scaffolds.csv` | ~14,000 |
| 5 | `scripts/phase1_data/fix3_acceptor_filter.py` | `data/interim/fix2_with_scaffolds.csv` | `data/interim/fix3_acceptors_only.csv` | ~11,000 |
| 6 | `scripts/phase1_data/fix4_bootstrap_leakage.py` | `data/interim/fix3_acceptors_only.csv` | `data/processed/master_acceptor_dataset.csv` | ~11,000 |

## Phase 2: Causal Inference

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 7 | `scripts/phase2_causal/dml_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/dml_ewg_corrected.csv` |
| 8 | `scripts/phase2_causal/causal_forest.R` | `data/processed/master_acceptor_dataset.csv` | `data/processed/experimental_with_cates.csv`, `results/tables/cate_results_corrected.csv` |
| 9 | `scripts/phase2_causal/iv_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/iv_results.csv` |
| 10 | `scripts/phase2_causal/sensitivity_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/sensitivity_results.csv`, `results/figures/sensitivity_*.png` |

## Phase 3: Counterfactual Design

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 11 | `scripts/phase3_design/step1_counterfactual_design.py` | `data/processed/experimental_with_cates.csv` | `results/tables/phase3_top10_candidates.csv` |
| 12 | `scripts/phase3_design/step2_xtb_validation.py` | `results/tables/phase3_top10_candidates.csv` | `results/tables/phase3_xtb_validation_final.csv` |

## Raw data notes

- `data/raw/Data_Merged_with_SMILES.xlsx`: 1,573 experimental acceptor molecules
- `data/raw/cep_sample_15000.csv`: 15,000-row reproducible sample from Harvard CEPDB, sampled with `random_state=42` / `set.seed(42)`

## Known issues

- `measurement_type` column values are `"experiment"` and `"DFT_B3LYP"` (not `"experimental"`). All Phase 2 R scripts filter on `"experiment"`.
- `iv_evaluation_table4.csv` has been removed. The IV results table is `results/tables/iv_results.csv`, produced by `iv_analysis.R` with `nonarom_cc` as the sole instrument. `halogen_count` is not used as an instrument (exclusion restriction violation — it is a component of `ewg_weighted`).
- xTB must be installed separately and available on PATH for `step2_xtb_validation.py` to run.
- DML cross-fitting uses scaffold-grouped folds (not random folds) so cross-fitting residuals are honest at scaffold-family level.