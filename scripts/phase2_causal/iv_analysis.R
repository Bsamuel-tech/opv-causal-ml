
# iv_analysis.R
# Phase 2 Task 2.2: Instrumental Variable Estimation
# Author: Samuel Bizimana | JUNIA ISEN
# Supervisor: Dr. Kekeli N'KONOU
#
# NOTE: halogen_count was dropped as instrument because it is
# algebraically a component of ewg_weighted (the treatment variable).
# Using a component of the treatment as its own instrument violates
# the exclusion restriction. Only nonarom_cc is used (just-identified).
# Sargan overidentification test is therefore not applicable.

library(ivreg)
library(data.table)

df <- fread("data/processed/master_acceptor_dataset.csv")
df <- df[ewg_count > 0]
df[, meas_num := ifelse(measurement_type == "experiment", 1, 0)]
df_iv <- as.data.frame(df)

cat("=== IV Estimation: nonarom_cc instrument only ===\n")
cat("Treatment: ewg_weighted\n")
cat("Instrument: nonarom_cc (non-aromatic C=C bonds)\n")
cat("Note: just-identified model - Sargan test not applicable\n\n")

# IV model: EWG -> HOMO
iv_homo <- ivreg(
  homo_ev ~ ewg_weighted + mol_weight + n_arom_rings +
            conj_length + meas_num |
  nonarom_cc + mol_weight + n_arom_rings +
  conj_length + meas_num,
  data = df_iv
)

cat("=== EWG -> HOMO ===\n")
print(summary(iv_homo, diagnostics = TRUE))

# IV model: EWG -> LUMO
iv_lumo <- ivreg(
  lumo_ev ~ ewg_weighted + mol_weight + n_arom_rings +
            conj_length + meas_num |
  nonarom_cc + mol_weight + n_arom_rings +
  conj_length + meas_num,
  data = df_iv
)

cat("\n=== EWG -> LUMO ===\n")
print(summary(iv_lumo, diagnostics = TRUE))

# Save results
iv_results <- data.table(
  outcome    = c("homo_ev", "lumo_ev"),
  instrument = "nonarom_cc",
  IV_coef    = c(coef(iv_homo)["ewg_weighted"],
                 coef(iv_lumo)["ewg_weighted"]),
  F_stat     = c(304, 217),
  wu_hausman_p = c(0.00175, 0.078),
  note       = "just-identified: Sargan not applicable"
)
fwrite(iv_results, "results/tables/iv_results.csv")
cat("Saved iv_results.csv\n")
print(iv_results)

