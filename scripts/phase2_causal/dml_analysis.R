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
#
# n_rep = 20 independent permutations are used so results are not
# dependent on a single fold draw. DoubleML aggregates by median.
N_REP <- 20

build_scaffold_split <- function(df, seed) {
  set.seed(seed)
  unique_scaffolds <- unique(df$scaffold)
  n_unique         <- length(unique_scaffolds)
  perm             <- sample(n_unique)
  fold_labels      <- ((seq_len(n_unique) - 1L) %% 5L) + 1L
  scaffold_fold    <- data.table(
    scaffold = unique_scaffolds[perm],
    fold_id  = fold_labels
  )
  df_tmp <- merge(df, scaffold_fold, by = "scaffold", all.x = TRUE)
  
  # Verify scaffold integrity
  check <- df_tmp[, .(n_folds_per_scaffold = uniqueN(fold_id)), by = scaffold]
  if (max(check$n_folds_per_scaffold) > 1)
    stop(paste("Scaffold spans multiple folds at seed", seed))
  
  fold_assignment <- df_tmp$fold_id
  if (any(table(fold_assignment) == 0))
    stop(paste("Empty fold detected at seed", seed))
  
  list(
    list(
      train_ids = lapply(1:5, function(k) which(fold_assignment != k)),
      test_ids  = lapply(1:5, function(k) which(fold_assignment == k))
    )
  )
}

# Build all N_REP splits upfront (seeds 1..N_REP)
all_splits <- lapply(seq_len(N_REP), function(i) build_scaffold_split(df, seed = i))

# Save fold diagnostic for seed 1 as representative
set.seed(1)
unique_scaffolds <- unique(df$scaffold)
n_unique         <- length(unique_scaffolds)
perm             <- sample(n_unique)
fold_labels      <- ((seq_len(n_unique) - 1L) %% 5L) + 1L
scaffold_fold_s1 <- data.table(scaffold = unique_scaffolds[perm], fold_id = fold_labels)
df_s1            <- merge(df, scaffold_fold_s1, by = "scaffold", all.x = TRUE)
diag_tbl         <- df_s1[, .(n_molecules = .N, n_scaffolds = uniqueN(scaffold)), by = fold_id]
setorder(diag_tbl, fold_id)
cat("\nFold diagnostic table (seed 1, representative):\n")
print(diag_tbl)
fwrite(diag_tbl, file.path(results_tables, "dml_fold_diagnostics.csv"))

# Train/test overlap guard (fold 1, seed 1)
fold_s1   <- df_s1$fold_id
test1     <- which(fold_s1 == 1)
train1    <- which(fold_s1 != 1)
if (length(intersect(test1, train1)) > 0) stop("Train/test overlap in fold 1 (seed 1)")
cat("Train/test overlap check fold 1 seed 1: OK\n")

# -- Confounders -----------------------------------------------------------
confounders <- c("mol_weight", "n_arom_rings", "conj_length", "halogen_count")
# measurement_type is constant in causal subset (all "experiment") -- excluded

# -- Helper to fit one DML model with n_rep splits ------------------------
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
  cat(sprintf("\n=== EWG -> %s (scaffold-grouped, %d reps) ===\n", label, N_REP))
  cat("ATE:", dml_obj$coef, "  SE:", dml_obj$se, "  p:", dml_obj$pval, "\n")
  
  # Extract per-rep ATEs to report spread
  rep_coefs <- as.numeric(dml_obj$all_coef)
  cat(sprintf("Per-rep ATE range: [%.4f, %.4f]  median: %.4f\n",
              min(rep_coefs), max(rep_coefs), median(rep_coefs)))
  dml_obj
}

dml_homo <- fit_dml("homo_ev",    "HOMO")
dml_lumo <- fit_dml("lumo_ev",    "LUMO")
dml_gap  <- fit_dml("bandgap_ev", "Bandgap")

# -- Save results ----------------------------------------------------------
# Extract per-rep spread
homo_reps <- as.numeric(dml_homo$all_coef)
lumo_reps <- as.numeric(dml_lumo$all_coef)
gap_reps  <- as.numeric(dml_gap$all_coef)

results <- data.table(
  treatment      = "ewg_weighted",
  outcome        = c("homo_ev", "lumo_ev", "bandgap_ev"),
  ATE            = c(dml_homo$coef,   dml_lumo$coef,   dml_gap$coef),
  SE             = c(dml_homo$se,     dml_lumo$se,     dml_gap$se),
  p_value        = c(dml_homo$pval,   dml_lumo$pval,   dml_gap$pval),
  CI_lower       = c(dml_homo$coef - 1.96*dml_homo$se,
                     dml_lumo$coef - 1.96*dml_lumo$se,
                     dml_gap$coef  - 1.96*dml_gap$se),
  CI_upper       = c(dml_homo$coef + 1.96*dml_homo$se,
                     dml_lumo$coef + 1.96*dml_lumo$se,
                     dml_gap$coef  + 1.96*dml_gap$se),
  rep_ATE_min    = c(min(homo_reps),  min(lumo_reps),  min(gap_reps)),
  rep_ATE_max    = c(max(homo_reps),  max(lumo_reps),  max(gap_reps)),
  rep_ATE_median = c(median(homo_reps), median(lumo_reps), median(gap_reps)),
  fold_scheme    = "scaffold-grouped-cyclic",
  n_rep          = N_REP,
  n_obs          = n_obs,
  n_scaffolds    = n_scaffolds,
  n_folds        = 5L
)

out_path <- file.path(results_tables, "dml_ewg_corrected.csv")
fwrite(results, out_path)
cat("\nSaved", out_path, "\n")
print(results)