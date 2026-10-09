library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)
library(here)

# -- Paths ------------------------------------------------------------------
data_processed <- here::here("data", "processed")
results_tables <- here::here("results", "tables")
dir.create(results_tables, recursive = TRUE, showWarnings = FALSE)

df <- fread(file.path(data_processed, "master_acceptor_dataset.csv"))

# Filter to causal estimation subset: experimental molecules with EWG groups
# measurement_type == "experiment" (not "experimental") -- verified from master dataset
df <- df[measurement_type == "experiment" & ewg_count > 0]
cat("Rows in causal estimation subset:", nrow(df), "\n")
cat("Unique scaffolds:", df[, uniqueN(scaffold)], "\n")

# Guard: abort if subset is empty or has missing scaffold values
if (nrow(df) == 0) stop("No rows match causal estimation filter -- check measurement_type values")
if (any(is.na(df$scaffold))) stop("NA scaffolds detected -- rerun fix2_scaffold_split.py")

n_obs       <- nrow(df)
n_scaffolds <- df[, uniqueN(scaffold)]

# -- Scaffold-grouped 5-fold assignment ------------------------------------
# Each scaffold is assigned to exactly ONE fold so structurally similar
# molecules never cross fold boundaries during cross-fitting.
#
# We use a deterministic round-robin (cyclic) assignment after a random
# permutation. This guarantees equal-ish fold sizes and no empty folds,
# unlike sample(1:5, replace=TRUE) which can in theory leave a fold empty.
set.seed(42)
unique_scaffolds <- unique(df$scaffold)
n_unique <- length(unique_scaffolds)
perm <- sample(n_unique)                               # random permutation, no replacement
fold_labels <- ((seq_len(n_unique) - 1L) %% 5L) + 1L  # cyclic: 1,2,3,4,5,1,2,3,...
scaffold_fold <- data.table(
  scaffold = unique_scaffolds[perm],
  fold_id  = fold_labels
)
df <- merge(df, scaffold_fold, by = "scaffold", all.x = TRUE)

# Verify scaffold integrity
check <- df[, .(n_folds_per_scaffold = uniqueN(fold_id)), by = scaffold]
max_folds_per_scaffold <- max(check$n_folds_per_scaffold)
cat("Max folds per scaffold (must be 1):", max_folds_per_scaffold, "\n")
if (max_folds_per_scaffold > 1) stop("Scaffold appears in more than one fold")

fold_assignment <- df$fold_id
cat("Fold distribution (scaffold-grouped, cyclic assignment):\n")
print(table(fold_assignment))

# Guard: no empty fold
fold_counts <- table(fold_assignment)
if (any(fold_counts == 0)) stop("Empty fold detected")

# Save diagnostic table (Arindam requirement: document fold distribution)
diag_tbl <- df[, .(n_molecules = .N, n_scaffolds = uniqueN(scaffold)), by = fold_id]
setorder(diag_tbl, fold_id)
cat("\nFold diagnostic table:\n")
print(diag_tbl)
fwrite(diag_tbl, file.path(results_tables, "dml_fold_diagnostics.csv"))

# Verify: test observations in fold k are not in fold k's training set
test1  <- which(fold_assignment == 1)
train1 <- which(fold_assignment != 1)
if (length(intersect(test1, train1)) > 0) stop("Train/test overlap in fold 1")
cat("Train/test overlap check fold 1: OK\n")

# -- Confounders -----------------------------------------------------------
confounders <- c("mol_weight", "n_arom_rings", "conj_length", "halogen_count")
# measurement_type is constant in causal subset (all "experiment") -- excluded

# Build the nested list DoubleML expects for set_sample_splitting()
custom_split <- list(
  list(
    train_ids = lapply(1:5, function(k) which(fold_assignment != k)),
    test_ids  = lapply(1:5, function(k) which(fold_assignment == k))
  )
)

# Guard: every fold's test set is non-empty
for (k in 1:5) {
  if (length(custom_split[[1]]$test_ids[[k]]) == 0)
    stop(paste("Empty test set for fold", k))
}

# -- Helper to fit one DML model ------------------------------------------
fit_dml <- function(y_col, label) {
  dml_data <- DoubleMLData$new(
    df, y_col = y_col, d_cols = "ewg_weighted", x_cols = confounders
  )
  set.seed(42)
  dml_obj <- DoubleMLPLR$new(
    dml_data,
    ml_l = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
    ml_m = lrn("regr.ranger", num.trees = 500, min.node.size = 5),
    n_folds               = 5,
    draw_sample_splitting = FALSE
  )
  dml_obj$set_sample_splitting(custom_split)
  dml_obj$fit()
  cat(sprintf("\n=== EWG -> %s (scaffold-grouped folds) ===\n", label))
  cat("ATE:", dml_obj$coef, "  SE:", dml_obj$se, "  p:", dml_obj$pval, "\n")
  dml_obj
}

dml_homo <- fit_dml("homo_ev",    "HOMO")
dml_lumo <- fit_dml("lumo_ev",    "LUMO")
dml_gap  <- fit_dml("bandgap_ev", "Bandgap")

# -- Save results (with n_obs, n_scaffolds, n_folds as required) ----------
results <- data.table(
  treatment   = "ewg_weighted",
  outcome     = c("homo_ev", "lumo_ev", "bandgap_ev"),
  ATE         = c(dml_homo$coef, dml_lumo$coef, dml_gap$coef),
  SE          = c(dml_homo$se,   dml_lumo$se,   dml_gap$se),
  p_value     = c(dml_homo$pval, dml_lumo$pval, dml_gap$pval),
  CI_lower    = c(dml_homo$coef - 1.96*dml_homo$se,
                  dml_lumo$coef - 1.96*dml_lumo$se,
                  dml_gap$coef  - 1.96*dml_gap$se),
  CI_upper    = c(dml_homo$coef + 1.96*dml_homo$se,
                  dml_lumo$coef + 1.96*dml_lumo$se,
                  dml_gap$coef  + 1.96*dml_gap$se),
  fold_scheme = "scaffold-grouped",
  n_obs       = n_obs,
  n_scaffolds = n_scaffolds,
  n_folds     = 5L
)

out_path <- file.path(results_tables, "dml_ewg_corrected.csv")
fwrite(results, out_path)
cat("\nSaved", out_path, "\n")
print(results)