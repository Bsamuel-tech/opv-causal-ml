
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

# OLS model for sensemakr
ols_homo <- lm(homo_ev ~ ewg_weighted + mol_weight + n_arom_rings +
               conj_length + logp + halogen_count + meas_num,
               data = df_sens)

# Sensitivity analysis
sens <- sensemakr(
  model            = ols_homo,
  treatment        = "ewg_weighted",
  benchmark_covars = "mol_weight",
  kd               = 1:3
)

cat("=== Sensitivity Analysis Results ===\n")
summary(sens)

# Save contour plot
png("results/figures/sensitivity_ewg_homo.png",
    width = 800, height = 600, res = 120)
plot(sens)
dev.off()

# Save results
sens_results <- data.table(
  treatment        = "ewg_weighted",
  outcome          = "homo_ev",
  ols_coef         = coef(ols_homo)["ewg_weighted"],
  robustness_value = sens$sensitivity_stats$rv_q,
  rv_alpha05       = sens$sensitivity_stats$rv_qa,
  partial_r2       = sens$sensitivity_stats$r2yd_x,
  benchmark_covar  = "mol_weight"
)
fwrite(sens_results, "results/tables/sensitivity_results.csv")
cat("Saved sensitivity_results.csv\n")

