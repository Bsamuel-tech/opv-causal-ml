
# causal_forest.R
# Phase 2 Task 2.4: Causal Forest and CATE Maps
# Author: Samuel Bizimana | JUNIA ISEN
# Supervisor: Dr. Kekeli N'KONOU

library(grf)
library(data.table)
library(ggplot2)

df <- fread("data/processed/master_acceptor_dataset.csv")
df <- df[ewg_count > 0]
cat("Rows:", nrow(df), "\n")

# Feature matrix, treatment, outcome
X <- as.matrix(df[, .(mol_weight, n_arom_rings, conj_length, logp)])
W <- df$ewg_weighted
Y <- df$homo_ev

# Train causal forest
set.seed(42)
cat("Training causal forest (4000 trees)...\n")
cf <- causal_forest(
  X             = X,
  Y             = Y,
  W             = W,
  num.trees     = 4000,
  min.node.size = 5,
  seed          = 42
)

# Extract CATEs
cate_out <- predict(cf, estimate.variance = TRUE)
df[, CATE    := cate_out$predictions]
df[, CATE_se := sqrt(cate_out$variance.estimates)]

cat("\n=== CATE Summary ===\n")
cat("Mean CATE:", round(mean(df$CATE), 6), "\n")
cat("Min CATE: ", round(min(df$CATE), 6), "\n")
cat("Max CATE: ", round(max(df$CATE), 6), "\n")

# Calibration test
cal <- test_calibration(cf)
cat("\n=== Calibration Test ===\n")
print(cal)

# Save model
saveRDS(cf, "results/models/causal_forest_ewg_homo.rds")
cat("Saved causal_forest_ewg_homo.rds\n")

# Save CATE results
fwrite(df[, .(canonical_SMILES, homo_ev, ewg_weighted,
              mol_weight, CATE, CATE_se)],
       "results/tables/cate_results_corrected.csv")

# CATE distribution plot
p1 <- ggplot(df, aes(x = CATE)) +
  geom_histogram(bins = 60, fill = "#2166ac",
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept = 0, color = "red",
             linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = mean(df$CATE),
             color = "orange", linewidth = 1) +
  labs(
    title    = "Causal Effect of EWG Score on HOMO Energy",
    subtitle = "Conditional Average Treatment Effects",
    x        = "CATE (eV per unit Hammett-weighted EWG score)",
    y        = "Number of molecules",
    caption  = "Orange = mean CATE | Red = zero effect"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("results/figures/cate_ewg_homo_distribution.png",
       p1, width = 9, height = 5, dpi = 300)

# CATE vs molecular weight
p2 <- ggplot(df, aes(x = mol_weight, y = CATE)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "gray50") +
  geom_ribbon(aes(ymin = CATE - 1.96*CATE_se,
                  ymax = CATE + 1.96*CATE_se),
              alpha = 0.2, fill = "#2166ac") +
  geom_point(alpha = 0.3, size = 1, color = "#2166ac") +
  geom_smooth(method = "loess", color = "#1F4E79",
              linewidth = 1.2, se = FALSE) +
  labs(
    title   = "CATE vs Molecular Weight",
    x       = "Molecular Weight (g/mol)",
    y       = "CATE (eV per unit EWG score)",
    caption = "Shaded = 95% CI | Curve = LOESS smoother"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("results/figures/cate_ewg_vs_molweight.png",
       p2, width = 9, height = 5, dpi = 300)

cat("All figures saved.\n")

