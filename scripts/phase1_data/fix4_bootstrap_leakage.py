"""
Fix 15: Bootstrap leakage gap analysis
=======================================
Replaces single-point generalisation gap estimate with
bootstrapped distribution across multiple scaffold splits.

Author: Samuel Bizimana | JUNIA ISEN
Supervisor: Dr. Kekeli N'KONOU
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
import random

print("=== Fix 15: Bootstrap leakage gap analysis ===\n")

df = pd.read_csv("data/processed/master_acceptor_dataset.csv")
features = ['ewg_weighted', 'mol_weight', 'n_arom_rings',
            'conj_length', 'logp', 'halogen_count']
target   = 'homo_ev'

X = df[features].values
y = df[target].values

# Random split R2 (single estimate - stable)
from sklearn.model_selection import train_test_split
X_tr, X_te, y_tr, y_te = train_test_split(
    X, y, test_size=0.2, random_state=42)
rf = RandomForestRegressor(n_estimators=200, random_state=42)
rf.fit(X_tr, y_tr)
r2_random = r2_score(y_te, rf.predict(X_te))
print(f"Random split R2: {r2_random:.4f}")

# Bootstrap scaffold splits
print("Running 20 bootstrap scaffold splits...")
scaffolds = df['scaffold'].values
unique_scaffolds = list(set(scaffolds))

gaps = []
r2_scaffolds = []

for seed in range(20):
    random.seed(seed)
    shuffled = unique_scaffolds.copy()
    random.shuffle(shuffled)

    n = len(shuffled)
    test_scaffolds = set(shuffled[int(0.8*n):])

    train_idx = [i for i, s in enumerate(scaffolds)
                 if s not in test_scaffolds]
    test_idx  = [i for i, s in enumerate(scaffolds)
                 if s in test_scaffolds]

    if len(test_idx) < 10:
        continue

    rf_s = RandomForestRegressor(n_estimators=200,
                                  random_state=seed)
    rf_s.fit(X[train_idx], y[train_idx])
    r2_s = r2_score(y[test_idx], rf_s.predict(X[test_idx]))
    gap  = r2_random - r2_s

    r2_scaffolds.append(r2_s)
    gaps.append(gap)

gaps = np.array(gaps)
r2_scaffolds = np.array(r2_scaffolds)

print(f"\n=== Bootstrap Results (20 seeds) ===")
print(f"Random split R2:         {r2_random:.4f}")
print(f"Scaffold R2 mean:        {r2_scaffolds.mean():.4f}")
print(f"Scaffold R2 std:         {r2_scaffolds.std():.4f}")
print(f"Scaffold R2 95% CI:      [{np.percentile(r2_scaffolds,2.5):.4f}, "
      f"{np.percentile(r2_scaffolds,97.5):.4f}]")
print(f"Gap mean:                {gaps.mean():.4f}")
print(f"Gap std:                 {gaps.std():.4f}")
print(f"Gap 95% CI:              [{np.percentile(gaps,2.5):.4f}, "
      f"{np.percentile(gaps,97.5):.4f}]")
print(f"Leakage threshold:       0.15")
print(f"All gaps below 0.15:     {(gaps < 0.15).all()}")

# Save
results = pd.DataFrame({
    'seed':        list(range(len(gaps))),
    'r2_scaffold': r2_scaffolds,
    'gap':         gaps
})
results.to_csv(
    "results/tables/bootstrap_leakage_gap.csv",
    index=False
)

summary = pd.DataFrame({
    'metric': [
        'random_split_r2',
        'scaffold_r2_mean', 'scaffold_r2_std',
        'scaffold_r2_ci_lower', 'scaffold_r2_ci_upper',
        'gap_mean', 'gap_std',
        'gap_ci_lower', 'gap_ci_upper',
        'leakage_threshold', 'all_below_threshold'
    ],
    'value': [
        round(r2_random, 4),
        round(r2_scaffolds.mean(), 4),
        round(r2_scaffolds.std(), 4),
        round(np.percentile(r2_scaffolds, 2.5), 4),
        round(np.percentile(r2_scaffolds, 97.5), 4),
        round(gaps.mean(), 4),
        round(gaps.std(), 4),
        round(np.percentile(gaps, 2.5), 4),
        round(np.percentile(gaps, 97.5), 4),
        0.15,
        str((gaps < 0.15).all())
    ]
})
summary.to_csv(
    "results/tables/bootstrap_leakage_summary.csv",
    index=False
)
print("\nSaved bootstrap_leakage_gap.csv")
print("Saved bootstrap_leakage_summary.csv")