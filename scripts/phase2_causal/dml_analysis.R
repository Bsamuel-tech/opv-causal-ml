
library(DoubleML); library(mlr3); library(mlr3learners); library(data.table)
df <- fread("data/processed/master_dataset.csv")
df[, HOMO_LUMO_gap := HOMO - LUMO]
df[, outcome := PCE]

