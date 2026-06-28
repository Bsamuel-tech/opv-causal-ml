import pandas as pd
from rdkit import Chem

print("Fixing HOPV15 unit conversion and rebuilding master dataset...")

# ── 1. YOUR DATA ───────────────────────────────────────────────────
own = pd.read_excel("Data_Merged_with_SMILES.xlsx")
own = own.rename(columns={
    'SMILES_don':  'SMILES',
    'HOMO_D':      'HOMO',
    'LUMO_D':      'LUMO',
    'Voc':         'Voc',
    'Jsc':         'Jsc',
    'FF':          'fill_factor',
    'PCE':         'PCE'
})
own = own[['SMILES','HOMO','LUMO','Voc','Jsc','fill_factor','PCE']].copy()
own['source_db'] = 'your_data'
print(f"Your data: {len(own)} rows | HOMO mean: {own['HOMO'].mean():.3f} eV")

# ── 2. HOPV15 — re-read raw and convert properly ───────────────────
hopv_raw = pd.read_csv("hopv_raw.csv")
hopv_raw = hopv_raw.rename(columns={
    'ids': 'SMILES',
    'y1':  'HOMO_hartree',
    'y2':  'LUMO_hartree',
    'y5':  'PCE',
    'y6':  'Voc',
    'y7':  'Jsc',
    'y8':  'fill_factor'
})
# Convert Hartree to eV
hopv_raw['HOMO'] = hopv_raw['HOMO_hartree'] * 27.2114
hopv_raw['LUMO'] = hopv_raw['LUMO_hartree'] * 27.2114
hopv_raw = hopv_raw[['SMILES','HOMO','LUMO','Voc','Jsc','fill_factor','PCE']].copy()
hopv_raw['source_db'] = 'HOPV15'
hopv_raw = hopv_raw.dropna(subset=['HOMO','LUMO'])
# Filter to chemically sensible OPV range
hopv_raw = hopv_raw[(hopv_raw['HOMO'] > -12) & (hopv_raw['HOMO'] < -2)]
print(f"HOPV15:    {len(hopv_raw)} rows | HOMO mean: {hopv_raw['HOMO'].mean():.3f} eV")

# ── 3. CEP ─────────────────────────────────────────────────────────
cep_full = pd.read_csv("moldata.csv")
cep = cep_full.sample(n=12000, random_state=42).copy()
cep = cep.rename(columns={
    'SMILES_str':   'SMILES',
    'e_homo_alpha': 'HOMO',
    'e_lumo_alpha': 'LUMO',
    'pce':          'PCE',
    'voc':          'Voc',
    'jsc':          'Jsc'
})
cep['fill_factor'] = None
cep = cep[['SMILES','HOMO','LUMO','Voc','Jsc','fill_factor','PCE']].copy()
cep['source_db'] = 'CEP'
print(f"CEP:       {len(cep)} rows | HOMO mean: {cep['HOMO'].mean():.3f} eV")

# ── 4. MERGE + CANONICALIZE + DEDUPLICATE ──────────────────────────
master = pd.concat([own, hopv_raw, cep], ignore_index=True)

def canonical(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

master['canonical_SMILES'] = master['SMILES'].apply(canonical)
master = master.dropna(subset=['canonical_SMILES'])
master = master.drop_duplicates(subset=['canonical_SMILES'])

# Final sanity filter — remove physically impossible values
master = master[(master['HOMO'] > -12) & (master['HOMO'] < -2)]
master = master[(master['LUMO'] > -8) & (master['LUMO'] < 2)]
master = master[(master['PCE'] >= 0) & (master['PCE'] <= 25)]

master.to_csv("master_dataset.csv", index=False)

print(f"\n✅ master_dataset.csv rebuilt")
print(f"   Total rows: {len(master)}")
print(f"\nSource breakdown:")
print(master['source_db'].value_counts())
print(f"\nHOMO range by source:")
print(master.groupby('source_db')['HOMO'].agg(['min','max','mean']).round(3))
print(f"\nPCE range by source:")
print(master.groupby('source_db')['PCE'].agg(['min','max','mean']).round(3))