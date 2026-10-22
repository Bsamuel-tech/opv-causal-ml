#!/usr/bin/env bash
# Re-generates MD5 checksums for all files listed in the PIPELINE.md manifest.
# Run from the repository root after any pipeline rerun.
# Usage: bash scripts/generate_manifest_md5.sh

set -euo pipefail

FILES=(
  "data/raw/Data_Merged_with_SMILES.xlsx"
  "data/interim/step1_experimental_raw.csv"
  "data/interim/fix3_acceptors_only.csv"
  "data/processed/master_acceptor_dataset.csv"
  "results/tables/dml_ewg_corrected.csv"
  "results/tables/dml_cluster_robust.csv"
  "results/tables/dml_fold_diagnostics.csv"
  "results/tables/iv_results.csv"
  "results/tables/iv_diagnostics.csv"
  "results/tables/cate_results_corrected.csv"
  "data/processed/experimental_with_cates.csv"
  "results/tables/bootstrap_leakage_summary.csv"
  "results/tables/bootstrap_leakage_gap.csv"
  "results/tables/phase3_xtb_validation_final.csv"
  "results/tables/supplementary_information.csv"
)

echo "File | MD5"
echo "-----|----"
for f in "${FILES[@]}"; do
  if [ -f "$f" ]; then
    hash=$(md5sum "$f" | awk '{print $1}')
    echo "$f | $hash"
  else
    echo "$f | MISSING"
  fi
done
