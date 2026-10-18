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

df <- df[measurement_type == "experiment" & ewg_count > 0]
cat("Rows in causal estimation subset:", nrow(df), "\n")
cat("Unique scaffolds:", df[, uniqueN(scaffold)], "\n")

if (nrow(df) == 0) stop("No rows match causal estimation filter")
if (any(is.na(df$scaffold))) stop("NA scaffolds detected")

n_obs       <- nrow(df)
n_scaffolds <- df[, uniqueN(scaffold)]

# -- Scaffold-grouped fold assignment (cyclic, named-vector lookup) --------
# Each scaffold is assigned to exactly one fold via a round-robin cyclic
# assignment after a random permutation. Named-vector lookup keeps
# fold_assignment aligned to df's original row order (no merge reordering).
# N_REP independent permutations are used; DoubleML aggregates by median.

N_REP <- 20

build_scaffold_split <- function(df, seed) {
  set.seed(seed)
  unique_scaffolds <- unique(df$scaffold)
  n_unique         <- length(unique_scaffolds)
  perm             <- sample(n_unique)
  fold_labels      <- ((seq_len(n_unique) - 1L) %% 5L) + 1L
  
  scaffold_to_fold <- setNames(fold_labels, unique_scaffolds[perm])
  fold_assignment  <- scaffold_to_fold[df$scaffold]
  
  check <- tapply(fold_assignment, df$scaffold, function(x) length(unique(x)))
  if (max(check) > 1)
    stop(paste("Scaffold spans multiple folds at seed", seed))
  if (any(table(fold_assignment) == 0))
    stop(paste("Empty fold detected at seed", seed))
  
  list(
    train_ids = lapply(1:5, function(k) which(fold_assignment != k)),
    test_ids  = lapply(1:5, function(k) which(fold_assignment == k))
  )
}

all_splits <- lapply(seq_len(N_REP), function(i) build_scaffold_split(df, seed = i))

# -- Fold diagnostic (seed 1, representative) ------------------------------
seed1_folds <- build_scaffold_split(df, seed = 1)

fold_s1 <- integer(nrow(df))
for (k in 1:5) fold_s1[ seed1_folds$test_ids[[k]] ] <- k

diag_tbl <- data.table(scaffold = df$scaffold, fold_id = fold_s1)
diag_tbl <- diag_tbl[, .(n_molecules = .N, n_scaffolds = uniqueN(scaffold)), by = fold_id]
setorder(diag_tbl, fold_id)
cat("\nFold diagnostic (seed 1):\n")
print(diag_tbl)
fwrite(diag_tbl, file.path(results_tables, "dml_fold_diagnostics.csv"))

# Train/test overlap guard (fold 1, seed 1)
test1  <- seed1_folds$test_ids[[1]]
train1 <- seed1_folds$train_ids[[1]]
if (length(intersect(test1, train1)) > 0) stop("Train/test overlap in fold 1 seed 1")
cat("Overlap check: OK\n")

# -- Confounders -----------------------------------------------------------
# measurement_type is constant in causal subset (all "experiment") -- excluded
confounders <- c("mol_weight", "n_arom_rings", "conj_length", "halogen_count")

# -- Fit DML ---------------------------------------------------------------
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
    n_rep                 = N_REP,
    draw_sample_splitting = FALSE
  )
  dml_obj$set_sample_splitting(all_splits)
  dml_obj$fit()
  
  rep_coefs <- as.numeric(dml_obj$all_coef)
  cat(sprintf("\n=== EWG -> %s (%d reps) ===\n", label, N_REP))
  cat("Median ATE:", median(rep_coefs), "\n")
  cat("SE:", dml_obj$se, "  p:", dml_obj$pval, "\n")
  cat(sprintf("Rep range [p5, p95]: [%.4f, %.4f]\n",
              quantile(rep_coefs, 0.05), quantile(rep_coefs, 0.95)))
  dml_obj
}

dml_homo <- fit_dml("homo_ev",    "HOMO")
dml_lumo <- fit_dml("lumo_ev",    "LUMO")
dml_gap  <- fit_dml("bandgap_ev", "Bandgap")

# -- Results table ---------------------------------------------------------
homo_reps <- as.numeric(dml_homo$all_coef)
lumo_reps <- as.numeric(dml_lumo$all_coef)
gap_reps  <- as.numeric(dml_gap$all_coef)

results <- data.table(
  treatment      = "ewg_weighted",
  outcome        = c("homo_ev", "lumo_ev", "bandgap_ev"),
  ATE            = c(dml_homo$coef,   dml_lumo$coef,   dml_gap$coef),
  SE             = c(dml_homo$se,     dml_lumo$se,     dml_gap$se),
  p_value        = c(dml_homo$pval,   dml_lumo$pval,   dml_gap$pval),
  CI_lower       = c(dml_homo$coef - 1.96 * dml_homo$se,
                     dml_lumo$coef - 1.96 * dml_lumo$se,
                     dml_gap$coef  - 1.96 * dml_gap$se),
  CI_upper       = c(dml_homo$coef + 1.96 * dml_homo$se,
                     dml_lumo$coef + 1.96 * dml_lumo$se,
                     dml_gap$coef  + 1.96 * dml_gap$se),
  rep_ATE_p5     = c(quantile(homo_reps, 0.05),
                     quantile(lumo_reps, 0.05),
                     quantile(gap_reps,  0.05)),
  rep_ATE_median = c(median(homo_reps), median(lumo_reps), median(gap_reps)),
  rep_ATE_p95    = c(quantile(homo_reps, 0.95),
                     quantile(lumo_reps, 0.95),
                     quantile(gap_reps,  0.95)),
  fold_scheme    = "scaffold-grouped-cyclic",
  n_rep          = N_REP,
  n_obs          = n_obs,
  n_scaffolds    = n_scaffolds,
  n_folds        = 5L
)

out_path <- file.path(results_tables, "dml_ewg_corrected.csv")
fwrite(results, out_path)
cat("\nSaved:", out_path, "\n")
print(results)

# -- Post-hoc cluster-robust variance correction ---------------------------
# DoubleML PLR does not natively support clustered inference. This block
# applies a manual sandwich correction using the score function residuals
# from the final repetition, grouped by scaffold.
#
# Method: HC-cluster sandwich estimator on the PLR score.
# Reference: Cameron & Miller (2015), "A Practitioner's Guide to
# Cluster-Robust Inference", Journal of Human Resources 50(2), 317-372.

cluster_robust_se <- function(dml_obj, df, cluster_var = "scaffold") {
  theta  <- as.numeric(dml_obj$coef)
  psi_b  <- as.numeric(dml_obj$psi_b[, 1, N_REP])
  psi_a  <- as.numeric(dml_obj$psi_a[, 1, N_REP])
  
  clusters <- df[[cluster_var]]
  n        <- length(psi_b)
  
  score        <- psi_a * theta + psi_b
  cluster_ids  <- unique(clusters)
  G            <- length(cluster_ids)
  
  meat <- sum(sapply(cluster_ids, function(g) {
    idx <- which(clusters == g)
    sum(score[idx])^2
  }))
  
  bread  <- sum(psi_a)
  cr_var <- (G / (G - 1)) * (n / (n - 1)) * meat / (bread^2)
  sqrt(cr_var)
}

cr_se_homo <- cluster_robust_se(dml_homo, df)
cr_se_lumo <- cluster_robust_se(dml_lumo, df)
cr_se_gap  <- cluster_robust_se(dml_gap,  df)

cr_results <- data.table(
  outcome     = c("homo_ev", "lumo_ev", "bandgap_ev"),
  ATE         = c(dml_homo$coef, dml_lumo$coef, dml_gap$coef),
  SE_iid      = c(dml_homo$se,   dml_lumo$se,   dml_gap$se),
  SE_cluster  = c(cr_se_homo,    cr_se_lumo,    cr_se_gap),
  p_cluster   = c(
    2 * pnorm(-abs(dml_homo$coef / cr_se_homo)),
    2 * pnorm(-abs(dml_lumo$coef / cr_se_lumo)),
    2 * pnorm(-abs(dml_gap$coef  / cr_se_gap))
  ),
  CI_lower_cr = c(dml_homo$coef - 1.96 * cr_se_homo,
                  dml_lumo$coef - 1.96 * cr_se_lumo,
                  dml_gap$coef  - 1.96 * cr_se_gap),
  CI_upper_cr = c(dml_homo$coef + 1.96 * cr_se_homo,
                  dml_lumo$coef + 1.96 * cr_se_lumo,
                  dml_gap$coef  + 1.96 * cr_se_gap),
  n_clusters  = df[, uniqueN(scaffold)]
)

cat("\nCluster-robust SE correction (scaffold clusters):\n")
print(cr_results)

cr_path <- file.path(results_tables, "dml_cluster_robust.csv")
fwrite(cr_results, cr_path)
cat("Saved:", cr_path, "\n")