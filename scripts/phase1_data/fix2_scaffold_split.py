import pandas as pd
import numpy as np
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
import random

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT         = Path(__file__).resolve().parent.parent.parent
DATA_INTERIM = ROOT / "data" / "interim"
DATA_INTERIM.mkdir(parents=True, exist_ok=True)

print("=== Fix 2: Scaffold-based split ===\n")

df = pd.read_csv(DATA_INTERIM / "fix1_with_columns.csv")
print(f"Loaded: {len(df)} rows from fix1_with_columns.csv")

# ── 1. Compute Murcko scaffold for every molecule ─────────────────
def get_scaffold(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return None
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold)
    except:
        return None

print("Computing Murcko scaffolds...")
df['scaffold'] = df['canonical_SMILES'].apply(get_scaffold)
df = df.dropna(subset=['scaffold'])
print(f"Rows after scaffold computation: {len(df)}")
print(f"Unique scaffolds: {df['scaffold'].nunique()}")

# ── 2. Assign scaffolds to train/val/test ─────────────────────────
random.seed(42)
scaffolds = list(df['scaffold'].unique())
random.shuffle(scaffolds)

n = len(scaffolds)
train_end = int(0.70 * n)
val_end   = int(0.85 * n)

train_scaffolds = set(scaffolds[:train_end])
val_scaffolds   = set(scaffolds[train_end:val_end])
test_scaffolds  = set(scaffolds[val_end:])

def assign_split(scaffold):
    if scaffold in train_scaffolds:
        return 'train'
    elif scaffold in val_scaffolds:
        return 'val'
    else:
        return 'test'

df['split'] = df['scaffold'].apply(assign_split)

print(f"\nSplit breakdown:")
print(df['split'].value_counts())
print(f"\nVerification — scaffold overlap between train and test:")
train_scaf = set(df[df['split']=='train']['scaffold'])
test_scaf  = set(df[df['split']=='test']['scaffold'])
overlap = train_scaf.intersection(test_scaf)
print(f"Overlapping scaffolds: {len(overlap)} (must be 0)")

# ── 3. Compare random split vs scaffold split R² ──────────────────
features = ['ewg_weighted', 'mol_weight', 'n_arom_rings',
            'conj_length', 'logp', 'halogen_count']
target = 'homo_ev'

X = df[features].values
y = df[target].values

# Random split
X_tr, X_te, y_tr, y_te = train_test_split(X, y, test_size=0.2, random_state=42)
rf = RandomForestRegressor(n_estimators=200, random_state=42)
rf.fit(X_tr, y_tr)
r2_random = r2_score(y_te, rf.predict(X_te))

# Scaffold split
train_idx = df[df['split']=='train'].index
test_idx  = df[df['split']=='test'].index
X_tr_sc = df.loc[train_idx, features].values
y_tr_sc = df.loc[train_idx, target].values
X_te_sc = df.loc[test_idx, features].values
y_te_sc = df.loc[test_idx, target].values
rf2 = RandomForestRegressor(n_estimators=200, random_state=42)
rf2.fit(X_tr_sc, y_tr_sc)
r2_scaffold = r2_score(y_te_sc, rf2.predict(X_te_sc))

print(f"\n=== R² Comparison ===")
print(f"Random split R²:    {r2_random:.4f}")
print(f"Scaffold split R²:  {r2_scaffold:.4f}")
print(f"Generalisation gap: {r2_random - r2_scaffold:.4f}")

if r2_random - r2_scaffold > 0.15:
    print("Data leakage CONFIRMED — gap > 0.15")
else:
    print("No significant data leakage detected.")

# ── 4. Save interim snapshot ───────────────────────────────────────
out = DATA_INTERIM / "fix2_with_scaffolds.csv"
df.to_csv(out, index=False)
print(f"\nSaved: {out}")
