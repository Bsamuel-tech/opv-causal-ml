
library(ivreg)
library(data.table)
df <- fread("data/processed/master_acceptor_dataset.csv")
df <- df[ewg_count > 0]
df[, meas_num := ifelse(measurement_type == "experiment", 1, 0)]
df_iv <- as.data.frame(df)
iv_homo <- ivreg(
  homo_ev ~ ewg_weighted + mol_weight + n_arom_rings + conj_length + meas_num |
  halogen_count + nonarom_cc + mol_weight + n_arom_rings + conj_length + meas_num,
  data = df_iv)
print(summary(iv_homo, diagnostics = TRUE))

