
# sensitivity_analysis.R
# Phase 2 Task 2.3: Sensitivity Analysis using sensemakr
# Author: Samuel Bizimana | JUNIA ISEN
# Supervisor: Dr. Kekeli N'KONOU

library(sensemakr)
library(data.table)

df <- fread("data/processed/master_acceptor_dataset.csv")
df <- df[ewg_count > 0]
df[, meas_num := ifelse(measurement_type == "experiment", 1, 0)]
df_sens <- as.data.frame(df)

# === PRIMARY: OLS model for LUMO ===
ols_lumo <- lm(lumo_ev ~ ewg_weighted + mol_weight + n_arom_rings +
               conj_length + halogen_count + meas_num,
               data = df_sens)

sens_lumo <- sensemakr(
  model            = ols_lumo,
  treatment        = "ewg_weighted",
  benchmark_covars = "mol_weight",
  kd               = 1:3
)

cat("=== PRIMARY: Sensitivity Analysis — EWG -> LUMO ===
")
summary(sens_lumo)

png("results/figures/sensitivity_ewg_lumo.png",
    width = 800, height = 600, res = 120)
plot(sens_lumo,
     main = "Sensitivity Analysis: EWG Score -> LUMO Energy (Primary)")
dev.off()
cat("Saved sensitivity_ewg_lumo.png
")

# === SECONDARY: OLS model for HOMO ===
ols_homo <- lm(homo_ev ~ ewg_weighted + mol_weight + n_arom_rings +
               conj_length + halogen_count + meas_num,
               data = df_sens)

sens_homo <- sensemakr(
  model            = ols_homo,
  treatment        = "ewg_weighted",
  benchmark_covars = "mol_weight",
  kd               = 1:3
)

cat("
=== SECONDARY: Sensitivity Analysis — EWG -> HOMO ===
")
summary(sens_homo)

png("results/figures/sensitivity_ewg_homo.png",
    width = 800, height = 600, res = 120)
plot(sens_homo,
     main = "Sensitivity Analysis: EWG Score -> HOMO Energy (Secondary)")
dev.off()
cat("Saved sensitivity_ewg_homo.png
")

# === Save results ===
sens_results <- data.table(
  outcome          = c("lumo_ev", "homo_ev"),
  ols_coef         = c(coef(ols_lumo)["ewg_weighted"],
                       coef(ols_homo)["ewg_weighted"]),
  robustness_value = c(sens_lumo$sensitivity_stats$rv_q,
                       sens_homo$sensitivity_stats$rv_q),
  rv_alpha05       = c(sens_lumo$sensitivity_stats$rv_qa,
                       sens_homo$sensitivity_stats$rv_qa),
  benchmark_covar  = "mol_weight",
  priority         = c("PRIMARY", "SECONDARY")
)
fwrite(sens_results, "results/tables/sensitivity_results.csv")
cat("Saved sensitivity_results.csv
")
print(sens_results)

