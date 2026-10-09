import pandas as pd
from pathlib import Path
from rdkit import Chem

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT         = Path(__file__).resolve().parent.parent.parent
DATA_INTERIM = ROOT / "data" / "interim"
DATA_INTERIM.mkdir(parents=True, exist_ok=True)

print("=== Fix 1: Adding measurement_type, halogen_count, nonarom_cc ===\n")

df = pd.read_csv(DATA_INTERIM / "step2_merged_raw.csv")
print(f"Loaded: {len(df)} rows from step2_merged_raw.csv")

# ── measurement_type ───────────────────────────────────────────────
df['measurement_type'] = df['source_db'].map({
    'your_data': 'experiment',
    'CEP':       'DFT_B3LYP'
})
print("measurement_type added:")
print(df['measurement_type'].value_counts())

# ── halogen_count ──────────────────────────────────────────────────
def count_halogens(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return 0
        return len(mol.GetSubstructMatches(
            Chem.MolFromSmarts('[F,Cl,Br,I]')
        ))
    except:
        return 0

print("\nComputing halogen counts...")
df['halogen_count'] = df['canonical_SMILES'].apply(count_halogens)
print(f"Halogen count range: {df['halogen_count'].min()} to {df['halogen_count'].max()}")
print(f"Molecules with halogens: {(df['halogen_count'] > 0).sum()}")

# ── non-aromatic C=C bonds (IV instrument) ─────────────────────────
def count_nonarom_cc(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return 0
        pattern = Chem.MolFromSmarts('[C;!a]=[C;!a]')
        return len(mol.GetSubstructMatches(pattern))
    except:
        return 0

print("Computing non-aromatic C=C bonds...")
df['nonarom_cc'] = df['canonical_SMILES'].apply(count_nonarom_cc)
print(f"Non-arom C=C range: {df['nonarom_cc'].min()} to {df['nonarom_cc'].max()}")

# ── Save interim snapshot ──────────────────────────────────────────
out = DATA_INTERIM / "fix1_with_columns.csv"
df.to_csv(out, index=False)
print(f"\nSaved: {out}")
print(f"Columns now: {df.columns.tolist()}")
print(f"Total rows: {len(df)}")

