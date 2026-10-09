library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)
library(here)

# ── Paths ──────────────────────────────────────────────────────────────────
data_processed <- here::here("data", "processed")
results_tables <- here::here("results", "tables")
dir.create(results_tables, recursive = TRUE, showWarnings = FALSE)

df <- fread(file.path(data_processed, "master_acceptor_dataset.csv"))
df <- df[ewg_count > 0]
cat("Rows in analysis subset:", nrow(df), "\n")
cat("Unique scaffolds:", df[, uniqueN(scaffold)], "\n")

# ── Scaffold-grouped 5-fold assignment ────────────────────────────────────
# Assign each scaffold to one of 5 folds so structurally similar
# molecules never split across fold boundaries — connecting the
# Phase 1 scaffold split to the DML cross-fitting step.
set.seed(42)
unique_scaffolds <- unique(df$scaffold)
scaffold_fold    <- data.table(
  scaffold   = unique_scaffolds,
  fold_id    = sample(1:5, length(unique_scaffolds), replace = TRUE)
)
df <- merge(df, scaffold_fold, by = "scaffold", all.x = TRUE)

# Build the fold assignment vector DoubleML expects (integer 1..n_folds)
fold_assignment <- df$fold_id
cat("Fold distribution (scaffold-grouped):\n")
print(table(fold_assignment))

# Verify: no scaffold spans two folds
check <- df[, .(n_folds = uniqueN(fold_id)), by = scaffold]
cat("Max folds per scaffold (must be 1):", max(check$n_folds), "\n")

# ── Confounders ───────────────────────────────────────────────────────────
confounders <- c("mol_weight", "n_arom_rings", "conj_length", "halogen_count")
# Note: meas_num excluded — constant across causal subset (all ewg_count > 0)

# ── Model 1: EWG -> HOMO ─────────────────────────────────────────────────
dml_data_homo <- DoubleMLData$new(
  df, y_col = "homo_ev", d_cols = "ewg_weighted", x_cols = confounders
)
set.seed(42)
dml_homo <- DoubleMLPLR$new(
  dml_data_homo,
  ml_l = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
  ml_m = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
  n_folds    = 5,
  draw_sample_splitting = FALSE   # we supply our own split below
)
dml_homo$set_sample_splitting(list(list(train_ids = lapply(1:5, function(k) which(fold_assignment != k)),
                                        test_ids  = lapply(1:5, function(k) which(fold_assignment == k)))))
dml_homo$fit()
cat("\n=== EWG -> HOMO (scaffold-grouped folds) ===\n")
cat("ATE:", dml_homo$coef, "  SE:", dml_homo$se, "  p:", dml_homo$pval, "\n")

# ── Model 2: EWG -> LUMO ─────────────────────────────────────────────────
dml_data_lumo <- DoubleMLData$new(
  df, y_col = "lumo_ev", d_cols = "ewg_weighted", x_cols = confounders
)
set.seed(42)
dml_lumo <- DoubleMLPLR$new(
  dml_data_lumo,
  ml_l = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
  ml_m = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
  n_folds    = 5,
  draw_sample_splitting = FALSE
)
dml_lumo$set_sample_splitting(list(list(train_ids = lapply(1:5, function(k) which(fold_assignment != k)),
                                        test_ids  = lapply(1:5, function(k) which(fold_assignment == k)))))
dml_lumo$fit()
cat("\n=== EWG -> LUMO (scaffold-grouped folds) ===\n")
cat("ATE:", dml_lumo$coef, "  SE:", dml_lumo$se, "  p:", dml_lumo$pval, "\n")

# ── Model 3: EWG -> Bandgap ──────────────────────────────────────────────
dml_data_gap <- DoubleMLData$new(
  df, y_col = "bandgap_ev", d_cols = "ewg_weighted", x_cols = confounders
)
set.seed(42)
dml_gap <- DoubleMLPLR$new(
  dml_data_gap,
  ml_l = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
  ml_m = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
  n_folds    = 5,
  draw_sample_splitting = FALSE
)
dml_gap$set_sample_splitting(list(list(train_ids = lapply(1:5, function(k) which(fold_assignment != k)),
                                       test_ids  = lapply(1:5, function(k) which(fold_assignment == k)))))
dml_gap$fit()
cat("\n=== EWG -> Bandgap (scaffold-grouped folds) ===\n")
cat("ATE:", dml_gap$coef, "  SE:", dml_gap$se, "  p:", dml_gap$pval, "\n")

# ── Save results ──────────────────────────────────────────────────────────
results <- data.table(
  treatment    = "ewg_weighted",
  outcome      = c("homo_ev", "lumo_ev", "bandgap_ev"),
  ATE          = c(dml_homo$coef, dml_lumo$coef, dml_gap$coef),
  SE           = c(dml_homo$se,   dml_lumo$se,   dml_gap$se),
  p_value      = c(dml_homo$pval, dml_lumo$pval, dml_gap$pval),
  fold_scheme  = "scaffold-grouped"
)
out_path <- file.path(results_tables, "dml_ewg_corrected.csv")
fwrite(results, out_path)
cat("\nSaved", out_path, "\n")
print(results)