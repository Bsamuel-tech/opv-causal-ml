import pandas as pd
from rdkit import Chem

print("=== Building master dataset ===\n")

# ── 1. YOUR EXISTING DATA ──────────────────────────────────────────
print("Loading your existing OPV data...")
own = pd.read_excel("Data_Merged_with_SMILES.xlsx")
own = own.rename(columns={
    'SMILES_don': 'SMILES',
    'HOMO_D':     'HOMO',
    'LUMO_D':     'LUMO',
    'Voc':        'Voc',
    'Jsc':        'Jsc',
    'FF':         'fill_factor',
    'PCE':        'PCE'
})
own = own[['SMILES','HOMO','LUMO','Voc','Jsc','fill_factor','PCE']].copy()
own['source_db'] = 'your_data'
print(f"  Your data: {len(own)} rows")

# ── 2. HOPV15 ──────────────────────────────────────────────────────
print("Loading HOPV15...")
hopv = pd.read_csv("hopv_clean.csv")
# Convert Hartree to eV (1 Hartree = 27.2114 eV)
hopv['HOMO'] = hopv['HOMO'] * 27.2114
hopv['LUMO'] = hopv['LUMO'] * 27.2114
hopv = hopv[['SMILES','HOMO','LUMO','Voc','Jsc','fill_factor','PCE']].copy()
hopv['source_db'] = 'HOPV15'
print(f"  HOPV15: {len(hopv)} rows")

# ── 3. CEP (sample 12,000 rows) ────────────────────────────────────
print("Loading CEP dataset (sampling 12,000 rows)...")
cep_full = pd.read_csv("moldata.csv")
cep = cep_full.sample(n=12000, random_state=42).copy()
cep = cep.rename(columns={
    'SMILES_str':    'SMILES',
    'e_homo_alpha':  'HOMO',
    'e_lumo_alpha':  'LUMO',
    'pce':           'PCE',
    'voc':           'Voc',
    'jsc':           'Jsc'
})
cep['fill_factor'] = None
cep = cep[['SMILES','HOMO','LUMO','Voc','Jsc','fill_factor','PCE']].copy()
cep['source_db'] = 'CEP'
print(f"  CEP sample: {len(cep)} rows")

# ── 4. MERGE ───────────────────────────────────────────────────────
print("\nMerging all sources...")
master = pd.concat([own, hopv, cep], ignore_index=True)
print(f"  Combined rows before dedup: {len(master)}")

# ── 5. CANONICALIZE & DEDUPLICATE ─────────────────────────────────
print("Canonicalizing SMILES and removing duplicates...")
def canonical(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

master['canonical_SMILES'] = master['SMILES'].apply(canonical)
master = master.dropna(subset=['canonical_SMILES'])
master = master.drop_duplicates(subset=['canonical_SMILES'])
print(f"  Rows after dedup: {len(master)}")

# ── 6. SAVE ────────────────────────────────────────────────────────
master.to_csv("master_dataset.csv", index=False)
print(f"\n✅ Saved master_dataset.csv")
print(f"   Total rows:    {len(master)}")
print(f"   Total columns: {len(master.columns)}")
print(f"\nSource breakdown:")
print(master['source_db'].value_counts())
print(f"\nMissing values:")
print(master[['HOMO','LUMO','PCE']].isnull().sum())