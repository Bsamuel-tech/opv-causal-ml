# PIPELINE.md — OPV Causal ML Pipeline
Run scripts in the order listed. Each step reads from the previous step's output.

## Phase 1 — Data

| Step | Script | Input | Output | Rows (approx) |
|------|--------|-------|--------|---------------|
| 1 | `scripts/phase1_data/step1_build_dataset.py` | `data/raw/Data_Merged_with_SMILES.xlsx` | `data/interim/step1_experimental_raw.csv` | 1,573 |
| 2 | `scripts/phase1_data/step2_expand_dataset.py` | `data/interim/step1_experimental_raw.csv` + `data/raw/moldata.csv` | `data/interim/step2_merged_raw.csv` | ~14,000 |
| 3 | `scripts/phase1_data/fix3_acceptor_filter.py` | `data/interim/step2_merged_raw.csv` | `data/interim/fix3_acceptors_only.csv` | ~11,000 |
| 4 | `scripts/phase1_data/fix4_bootstrap_leakage.py` | `data/interim/fix3_acceptors_only.csv` | `data/processed/master_acceptor_dataset.csv` | ~11,000 |

## Phase 2 — Causal Inference

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 5 | `scripts/phase2_causal/causal_forest.R` | `data/processed/master_acceptor_dataset.csv` | `data/processed/experimental_with_cates.csv`, `results/tables/cate_results_corrected.csv`, `results/figures/cate_*.png` |
| 6 | `scripts/phase2_causal/iv_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/iv_results.csv` |
| 7 | `scripts/phase2_causal/sensitivity_analysis.R` | `data/processed/master_acceptor_dataset.csv` | `results/tables/sensitivity_results.csv`, `results/figures/sensitivity_*.png` |

## Phase 3 — Counterfactual Design

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 8 | `scripts/phase3_design/step1_counterfactual_design.py` | `data/processed/experimental_with_cates.csv` | `results/tables/phase3_top10_candidates.csv` |
| 9 | `scripts/phase3_design/step2_xtb_validation.py` | `results/tables/phase3_top10_candidates.csv` | `results/tables/phase3_xtb_validation_final.csv` |

## Raw data notes

- `data/raw/Data_Merged_with_SMILES.xlsx` — 1,573 experimental acceptor molecules (not committed; attach manually)
- `data/raw/moldata.csv` — Harvard CEPDB; downloaded [DATE]; SHA256: [CHECKSUM]; step2 samples 15,000 rows with `random_state=42`

## Known issues

- `moldata.csv` is too large for GitHub. Either commit `data/raw/cep_sample_15000.csv` (the exact 15k rows used) or document the checksum above so others can reproduce the sample.
- xTB must be installed separately and available on PATH for step2_xtb_validation.py to run.