# iv_analysis.R
# Phase 2 Task 2.2: Instrumental Variable Estimation

# NOTE: halogen_count was dropped as instrument because it is
# algebraically a component of ewg_weighted (the treatment variable).
# Using a component of the treatment as its own instrument violates
# the exclusion restriction. Only nonarom_cc is used (just-identified).
# Sargan overidentification test is therefore not applicable.

library(here)
library(ivreg)
library(data.table)

# ── Paths ──────────────────────────────────────────────────────────────────
data_processed <- here::here("data", "processed")
results_tables <- here::here("results", "tables")
dir.create(results_tables, recursive = TRUE, showWarnings = FALSE)

df <- fread(file.path(data_processed, "master_acceptor_dataset.csv"))
df <- df[measurement_type == "experimental" & ewg_count > 0]
# NOTE: measurement_type is constant in causal subset (all experimental),
#       so meas_num is excluded from IV covariates — it adds no information.
df_iv <- as.data.frame(df)

cat("=== IV Estimation: nonarom_cc instrument only ===\n")
cat("Rows in causal subset:", nrow(df_iv), "\n")
cat("Treatment: ewg_weighted\n")
cat("Instrument: nonarom_cc (non-aromatic C=C bonds)\n")
cat("Note: just-identified model - Sargan test not applicable\n\n")

# IV model: EWG -> HOMO
iv_homo <- ivreg(
  homo_ev ~ ewg_weighted + mol_weight + n_arom_rings + conj_length |
    nonarom_cc + mol_weight + n_arom_rings + conj_length,
  data = df_iv
)

cat("=== EWG -> HOMO ===\n")
print(summary(iv_homo, diagnostics = TRUE))

# IV model: EWG -> LUMO
iv_lumo <- ivreg(
  lumo_ev ~ ewg_weighted + mol_weight + n_arom_rings + conj_length |
    nonarom_cc + mol_weight + n_arom_rings + conj_length,
  data = df_iv
)

cat("\n=== EWG -> LUMO ===\n")
print(summary(iv_lumo, diagnostics = TRUE))

# Save results
iv_results <- data.table(
  outcome      = c("homo_ev", "lumo_ev"),
  instrument   = "nonarom_cc",
  IV_coef      = c(coef(iv_homo)["ewg_weighted"],
                   coef(iv_lumo)["ewg_weighted"]),
  note         = "just-identified: Sargan not applicable; run diagnostics=TRUE for F-stat and Wu-Hausman"
)

out_path <- file.path(results_tables, "iv_results.csv")
fwrite(iv_results, out_path)
cat("Saved", out_path, "\n")
print(iv_results)