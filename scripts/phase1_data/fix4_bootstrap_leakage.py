"""
Fix 15: Bootstrap leakage gap analysis
=======================================
Replaces single-point generalisation gap estimate with
bootstrapped distribution across multiple scaffold splits.

IMPORTANT SCIENTIFIC NOTE (Arindam audit, 2026-10-09):
The full-dataset bootstrap (Section A below) uses ~15,000 molecules
from fix3_acceptors_only.csv, which is dominated by CEP DFT molecules.
The DML causal estimation (Phase 2) uses only the 598 experimental
molecules with ewg_count > 0 — a structurally different subset.
The full-dataset result cannot serve as proof that the causal subset
has no scaffold leakage. Section B runs a separate bootstrap on the
actual causal subset and is the relevant check for Phase 2 validity.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
import random

# -- Paths ------------------------------------------------------------------
ROOT           = Path(__file__).resolve().parent.parent.parent
DATA_INTERIM   = ROOT / "data" / "interim"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS_TABLES = ROOT / "results" / "tables"
RESULTS_TABLES.mkdir(parents=True, exist_ok=True)

print("=== Fix 15: Bootstrap leakage gap analysis ===\n")

df = pd.read_csv(DATA_INTERIM / "fix3_acceptors_only.csv")
print(f"Full dataset loaded: {len(df)} rows")

# Validate required columns
required_cols = ['canonical_SMILES', 'homo_ev', 'scaffold',
                 'ewg_weighted', 'mol_weight', 'n_arom_rings',
                 'conj_length', 'logp', 'halogen_count',
                 'measurement_type', 'ewg_count']
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Validate no NaN in scaffold or homo_ev
assert df['scaffold'].isnull().sum() == 0, "NaN scaffolds detected -- rerun fix2"
assert df['homo_ev'].isnull().sum() == 0, "NaN homo_ev detected"

features = ['ewg_weighted', 'mol_weight', 'n_arom_rings',
            'conj_length', 'logp', 'halogen_count']
target   = 'homo_ev'

# ===========================================================================
# SECTION A: Full dataset bootstrap (includes CEP DFT molecules)
# This characterises generalisation gap for the prediction model on all data.
# It is NOT the relevant check for Phase 2 causal estimation.
# ===========================================================================
print("\n=== SECTION A: Full dataset (~15k molecules) ===")
print("(dominantly CEP DFT -- not the causal estimation sample)")

X = df[features].values
y = df[target].values

X_tr, X_te, y_tr, y_te = train_test_split(
    X, y, test_size=0.2, random_state=42)
rf = RandomForestRegressor(n_estimators=200, random_state=42)
rf.fit(X_tr, y_tr)
r2_random_full = r2_score(y_te, rf.predict(X_te))
print(f"Random split R2 (full): {r2_random_full:.4f}")

scaffolds_full = df['scaffold'].values
unique_scaffolds_full = list(set(scaffolds_full))
gaps_full, r2_scaffolds_full = [], []

print("Running 20 bootstrap scaffold splits on full dataset...")
for seed in range(20):
    random.seed(seed)
    shuffled = unique_scaffolds_full.copy()
    random.shuffle(shuffled)
    n = len(shuffled)
    test_scaffolds = set(shuffled[int(0.8*n):])
    train_idx = [i for i, s in enumerate(scaffolds_full) if s not in test_scaffolds]
    test_idx  = [i for i, s in enumerate(scaffolds_full) if s in test_scaffolds]
    if len(test_idx) < 10:
        continue
    rf_s = RandomForestRegressor(n_estimators=200, random_state=seed)
    rf_s.fit(X[train_idx], y[train_idx])
    r2_s = r2_score(y[test_idx], rf_s.predict(X[test_idx]))
    gaps_full.append(r2_random_full - r2_s)
    r2_scaffolds_full.append(r2_s)

gaps_full         = np.array(gaps_full)
r2_scaffolds_full = np.array(r2_scaffolds_full)
print(f"Full dataset gap mean: {gaps_full.mean():.4f}  std: {gaps_full.std():.4f}")
print(f"All gaps below 0.15:  {(gaps_full < 0.15).all()}")

# ===========================================================================
# SECTION B: Causal subset bootstrap (experiment + ewg_count > 0)
# This is the relevant check for Phase 2 DML and causal forest validity.
# ===========================================================================
print("\n=== SECTION B: Causal subset (experiment + ewg_count > 0) ===")
print("(the 598 molecules actually used in DML / causal forest)")

causal = df[(df['measurement_type'] == 'experiment') & (df['ewg_count'] > 0)].copy()
print(f"Causal subset: {len(causal)} rows, {causal['scaffold'].nunique()} unique scaffolds")

if len(causal) < 50:
    print("WARNING: causal subset too small for bootstrap -- skipping Section B")
    gaps_causal         = np.array([])
    r2_scaffolds_causal = np.array([])
    r2_random_causal    = float('nan')
else:
    X_c = causal[features].values
    y_c = causal[target].values
    scaffolds_causal        = causal['scaffold'].values
    unique_scaffolds_causal = list(set(scaffolds_causal))

    X_tr_c, X_te_c, y_tr_c, y_te_c = train_test_split(
        X_c, y_c, test_size=0.2, random_state=42)
    rf_c = RandomForestRegressor(n_estimators=200, random_state=42)
    rf_c.fit(X_tr_c, y_tr_c)
    r2_random_causal = r2_score(y_te_c, rf_c.predict(X_te_c))
    print(f"Random split R2 (causal subset): {r2_random_causal:.4f}")

    gaps_causal, r2_scaffolds_causal = [], []
    print("Running 20 bootstrap scaffold splits on causal subset...")
    for seed in range(20):
        random.seed(seed)
        shuffled = unique_scaffolds_causal.copy()
        random.shuffle(shuffled)
        n = len(shuffled)
        test_scaffolds  = set(shuffled[int(0.8*n):])
        train_idx_c = [i for i, s in enumerate(scaffolds_causal) if s not in test_scaffolds]
        test_idx_c  = [i for i, s in enumerate(scaffolds_causal) if s in test_scaffolds]
        if len(test_idx_c) < 5:
            continue
        rf_cs = RandomForestRegressor(n_estimators=200, random_state=seed)
        rf_cs.fit(X_c[train_idx_c], y_c[train_idx_c])
        r2_cs = r2_score(y_c[test_idx_c], rf_cs.predict(X_c[test_idx_c]))
        gaps_causal.append(r2_random_causal - r2_cs)
        r2_scaffolds_causal.append(r2_cs)

    gaps_causal         = np.array(gaps_causal)
    r2_scaffolds_causal = np.array(r2_scaffolds_causal)
    print(f"Causal subset gap mean: {gaps_causal.mean():.4f}  std: {gaps_causal.std():.4f}")
    print(f"All gaps below 0.15:    {(gaps_causal < 0.15).all()}")
    print()
    print("=== Causal Subset Bootstrap Results (20 seeds) ===")
    print(f"Random split R2:         {r2_random_causal:.4f}")
    print(f"Scaffold R2 mean:        {r2_scaffolds_causal.mean():.4f}")
    print(f"Scaffold R2 std:         {r2_scaffolds_causal.std():.4f}")
    print(f"Gap mean:                {gaps_causal.mean():.4f}")
    print(f"Gap 95% CI:              [{np.percentile(gaps_causal,2.5):.4f}, "
          f"{np.percentile(gaps_causal,97.5):.4f}]")
    print(f"All gaps below 0.15:     {(gaps_causal < 0.15).all()}")

# -- Save per-seed results --------------------------------------------------
gap_out = RESULTS_TABLES / "bootstrap_leakage_gap.csv"
results = pd.DataFrame({
    'seed':               list(range(len(gaps_full))),
    'r2_scaffold_full':   r2_scaffolds_full,
    'r2_random_full':     r2_random_full,
    'gap_full':           gaps_full,
    'r2_scaffold_causal': r2_scaffolds_causal if len(r2_scaffolds_causal) == len(gaps_full) else [float('nan')]*len(gaps_full),
    'r2_random_causal':   r2_random_causal,
    'gap_causal':         gaps_causal if len(gaps_causal) == len(gaps_full) else [float('nan')]*len(gaps_full),
})
results.to_csv(gap_out, index=False)

# -- Save summary -----------------------------------------------------------
summary_out = RESULTS_TABLES / "bootstrap_leakage_summary.csv"
rows = [
    ('population',             'full_dataset'),
    ('n_molecules_full',       str(len(df))),
    ('random_split_r2_full',   round(r2_random_full, 4)),
    ('scaffold_r2_mean_full',  round(r2_scaffolds_full.mean(), 4)),
    ('scaffold_r2_std_full',   round(r2_scaffolds_full.std(), 4)),
    ('gap_mean_full',          round(gaps_full.mean(), 4)),
    ('gap_std_full',           round(gaps_full.std(), 4)),
    ('all_below_0.15_full',    str((gaps_full < 0.15).all())),
    ('population_causal',      'experiment+ewg_count>0'),
    ('n_molecules_causal',     str(len(causal))),
    ('n_scaffolds_causal',     str(causal['scaffold'].nunique())),
    ('random_split_r2_causal', round(r2_random_causal, 4) if not np.isnan(r2_random_causal) else 'NA'),
    ('scaffold_r2_mean_causal',round(r2_scaffolds_causal.mean(), 4) if len(r2_scaffolds_causal) > 0 else 'NA'),
    ('gap_mean_causal',        round(gaps_causal.mean(), 4) if len(gaps_causal) > 0 else 'NA'),
    ('gap_std_causal',         round(gaps_causal.std(), 4) if len(gaps_causal) > 0 else 'NA'),
    ('all_below_0.15_causal',  str((gaps_causal < 0.15).all()) if len(gaps_causal) > 0 else 'NA'),
    ('leakage_threshold',      0.15),
    ('note', 'Section A: full dataset (informative); Section B: causal subset (relevant for DML validity)'),
]
summary = pd.DataFrame(rows, columns=['metric', 'value'])
summary.to_csv(summary_out, index=False)

# Also copy final master dataset to processed/ for downstream scripts
master_out = DATA_PROCESSED / "master_acceptor_dataset.csv"
df.to_csv(master_out, index=False)

print(f"\nSaved: {gap_out}")
print(f"Saved: {summary_out}")
print(f"Saved: {master_out}")
print("\nDone. Use Section B (causal subset) results when discussing DML leakage.")
