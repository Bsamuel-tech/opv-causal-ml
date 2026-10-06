# sensitivity_analysis.R
# Phase 2 Task 2.3: Sensitivity Analysis using sensemakr

library(here)
library(sensemakr)
library(data.table)

# ── Paths ──────────────────────────────────────────────────────────────────
data_processed <- here::here("data", "processed")
results_tables <- here::here("results", "tables")
results_figures <- here::here("results", "figures")
dir.create(results_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(results_figures, recursive = TRUE, showWarnings = FALSE)

df <- fread(file.path(data_processed, "master_acceptor_dataset.csv"))
df <- df[measurement_type == "experimental" & ewg_count > 0]
# measurement_type is constant in causal subset — meas_num excluded
df_sens <- as.data.frame(df)

cat("Rows in causal subset:", nrow(df_sens), "\n")

# ── OLS: EWG -> HOMO ───────────────────────────────────────────────────────
ols_homo <- lm(homo_ev ~ ewg_weighted + mol_weight + n_arom_rings +
                 conj_length + halogen_count,
               data = df_sens)

sens_homo <- sensemakr(
  model            = ols_homo,
  treatment        = "ewg_weighted",
  benchmark_covars = "mol_weight",
  kd               = 1:3
)

cat("\n=== Sensitivity Analysis — EWG -> HOMO ===\n")
summary(sens_homo)

png(file.path(results_figures, "sensitivity_ewg_homo.png"),
    width = 800, height = 600, res = 120)
plot(sens_homo,
     main = "Sensitivity Analysis: EWG Score -> HOMO Energy")
dev.off()
cat("Saved sensitivity_ewg_homo.png\n")

# ── OLS: EWG -> LUMO ───────────────────────────────────────────────────────
ols_lumo <- lm(lumo_ev ~ ewg_weighted + mol_weight + n_arom_rings +
                 conj_length + halogen_count,
               data = df_sens)

sens_lumo <- sensemakr(
  model            = ols_lumo,
  treatment        = "ewg_weighted",
  benchmark_covars = "mol_weight",
  kd               = 1:3
)

cat("\n=== Sensitivity Analysis — EWG -> LUMO ===\n")
summary(sens_lumo)

png(file.path(results_figures, "sensitivity_ewg_lumo.png"),
    width = 800, height = 600, res = 120)
plot(sens_lumo,
     main = "Sensitivity Analysis: EWG Score -> LUMO Energy")
dev.off()
cat("Saved sensitivity_ewg_lumo.png\n")

# ── Save combined results ──────────────────────────────────────────────────
sens_results <- data.table(
  outcome          = c("homo_ev", "lumo_ev"),
  ols_coef         = c(coef(ols_homo)["ewg_weighted"],
                       coef(ols_lumo)["ewg_weighted"]),
  robustness_value = c(sens_homo$sensitivity_stats$rv_q,
                       sens_lumo$sensitivity_stats$rv_q),
  rv_alpha05       = c(sens_homo$sensitivity_stats$rv_qa,
                       sens_lumo$sensitivity_stats$rv_qa),
  benchmark_covar  = "mol_weight"
)

out_path <- file.path(results_tables, "sensitivity_results.csv")
fwrite(sens_results, out_path)
cat("Saved", out_path, "\n")
print(sens_results)