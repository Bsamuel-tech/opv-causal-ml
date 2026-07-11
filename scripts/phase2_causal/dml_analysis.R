
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)

df <- fread("data/processed/master_acceptor_dataset.csv")
df <- df[ewg_count > 0]
df[, meas_num := ifelse(measurement_type == "experiment", 1, 0)]
cat("Rows:", nrow(df), "
")

confounders <- c("mol_weight", "n_arom_rings", "conj_length",
                 "halogen_count", "meas_num")

# Model 1: EWG -> HOMO
dml_data_homo <- DoubleMLData$new(df, y_col="homo_ev",
  d_cols="ewg_weighted", x_cols=confounders)
set.seed(42)
dml_homo <- DoubleMLPLR$new(dml_data_homo,
  ml_l=lrn("regr.ranger", num.trees=500, min.node.size=5),
  ml_m=lrn("regr.ranger", num.trees=500, min.node.size=5))
dml_homo$fit()
cat("=== EWG -> HOMO ===
")
cat("ATE:", dml_homo$coef, "SE:", dml_homo$se, "p:", dml_homo$pval, "
")

# Model 2: EWG -> LUMO
dml_data_lumo <- DoubleMLData$new(df, y_col="lumo_ev",
  d_cols="ewg_weighted", x_cols=confounders)
set.seed(42)
dml_lumo <- DoubleMLPLR$new(dml_data_lumo,
  ml_l=lrn("regr.ranger", num.trees=500, min.node.size=5),
  ml_m=lrn("regr.ranger", num.trees=500, min.node.size=5))
dml_lumo$fit()
cat("=== EWG -> LUMO ===
")
cat("ATE:", dml_lumo$coef, "SE:", dml_lumo$se, "p:", dml_lumo$pval, "
")

# Model 3: EWG -> Bandgap
dml_data_gap <- DoubleMLData$new(df, y_col="bandgap_ev",
  d_cols="ewg_weighted", x_cols=confounders)
set.seed(42)
dml_gap <- DoubleMLPLR$new(dml_data_gap,
  ml_l=lrn("regr.ranger", num.trees=500, min.node.size=5),
  ml_m=lrn("regr.ranger", num.trees=500, min.node.size=5))
dml_gap$fit()
cat("=== EWG -> Bandgap ===
")
cat("ATE:", dml_gap$coef, "SE:", dml_gap$se, "p:", dml_gap$pval, "
")

# Save results
results <- data.table(
  treatment = "ewg_weighted",
  outcome   = c("homo_ev", "lumo_ev", "bandgap_ev"),
  ATE       = c(dml_homo$coef, dml_lumo$coef, dml_gap$coef),
  SE        = c(dml_homo$se,   dml_lumo$se,   dml_gap$se),
  p_value   = c(dml_homo$pval, dml_lumo$pval, dml_gap$pval)
)
fwrite(results, "results/tables/dml_ewg_corrected.csv")
cat("
Saved results/tables/dml_ewg_corrected.csv
")
print(results)

