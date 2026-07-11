import pandas as pd
from rdkit import Chem

print("=== Fix 1: Adding measurement_type and halogen_count ===\n")

df = pd.read_csv("master_acceptor_dataset.csv")
print(f"Loaded: {len(df)} rows")

# ── measurement_type ───────────────────────────────────────────────
# Your 599 experimental molecules vs CEP DFT molecules
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

# ── non-aromatic C=C bonds (second instrument for IV) ─────────────
def count_nonarom_cc(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return 0
        # Non-aromatic C=C double bonds
        pattern = Chem.MolFromSmarts('[C;!a]=[C;!a]')
        return len(mol.GetSubstructMatches(pattern))
    except:
        return 0

print("Computing non-aromatic C=C bonds...")
df['nonarom_cc'] = df['canonical_SMILES'].apply(count_nonarom_cc)
print(f"Non-arom C=C range: {df['nonarom_cc'].min()} to {df['nonarom_cc'].max()}")

# ── Save ───────────────────────────────────────────────────────────
df.to_csv("master_acceptor_dataset.csv", index=False)
print(f"\n✅ Saved master_acceptor_dataset.csv")
print(f"Columns now: {df.columns.tolist()}")
print(f"Total rows: {len(df)}")