import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

print("=== Step 2: Expanding to 10,000-15,000 acceptor molecules ===\n")

# ── 1. Load your base dataset ─────────────────────────────────────
base = pd.read_csv("acceptor_base.csv")
print(f"Base dataset: {len(base)} molecules")

# ── 2. Load CEP dataset ───────────────────────────────────────────
cep = pd.read_csv(
    r"C:\Users\Samuel Bizimana\OneDrive\Desktop\Research Training\data\raw\moldata.csv"
)
print(f"CEP raw: {len(cep)} rows")
print("CEP columns:", cep.columns.tolist())

# ── 3. Extract relevant CEP columns ───────────────────────────────
# CEP has: SMILES_str, e_homo_alpha, e_lumo_alpha, e_gap_alpha
cep = cep.rename(columns={
    'SMILES_str':    'SMILES',
    'e_homo_alpha':  'homo_ev',
    'e_lumo_alpha':  'lumo_ev',
    'e_gap_alpha':   'bandgap_ev'
})

# CEP energies are in eV already
cep = cep[['SMILES', 'homo_ev', 'lumo_ev', 'bandgap_ev']].copy()
cep['source_db'] = 'CEP'

# ── 4. Sample and clean CEP ───────────────────────────────────────
# Sample 15,000 to give room after filtering
cep_sample = cep.sample(n=15000, random_state=42).copy()

def canonical(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

print("Canonicalizing CEP SMILES...")
cep_sample['canonical_SMILES'] = cep_sample['SMILES'].apply(canonical)
cep_sample = cep_sample.dropna(subset=['canonical_SMILES'])
cep_sample = cep_sample.drop_duplicates(subset=['canonical_SMILES'])

# Physical range filter — same as your base data
cep_sample = cep_sample[(cep_sample['homo_ev'] > -8) & (cep_sample['homo_ev'] < -3)]
cep_sample = cep_sample[(cep_sample['lumo_ev'] > -6) & (cep_sample['lumo_ev'] < -1)]
cep_sample = cep_sample[(cep_sample['bandgap_ev'] > 0.5) & (cep_sample['bandgap_ev'] < 4.0)]
print(f"CEP after cleaning: {len(cep_sample)} molecules")

# ── 5. Compute EWG features for CEP ──────────────────────────────
def compute_ewg_features(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None, None, None, None, None, None
        cyano    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]#[N]')))
        carbonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=[O]')))
        nitro    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[N+](=O)[O-]')))
        halogen  = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F,Cl,Br,I]')))
        ester    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2][CX3](=[OX1])')))
        sulfonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[SX4](=[OX1])(=[OX1])')))
        ewg_count = cyano + carbonyl + nitro + halogen + ester + sulfonyl
        hammett = {
            'cyano': 0.66, 'carbonyl': 0.50, 'nitro': 0.78,
            'halogen': 0.23, 'ester': 0.45, 'sulfonyl': 0.72
        }
        ewg_weighted = (
            cyano    * hammett['cyano']    +
            carbonyl * hammett['carbonyl'] +
            nitro    * hammett['nitro']    +
            halogen  * hammett['halogen']  +
            ester    * hammett['ester']    +
            sulfonyl * hammett['sulfonyl']
        )
        mol_weight   = Descriptors.MolWt(mol)
        n_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        conj_length  = rdMolDescriptors.CalcNumRotatableBonds(mol)
        logp         = Descriptors.MolLogP(mol)
        return ewg_count, ewg_weighted, mol_weight, n_arom_rings, conj_length, logp
    except:
        return None, None, None, None, None, None

print("Computing EWG features for CEP molecules...")
results = cep_sample['canonical_SMILES'].apply(compute_ewg_features)
cep_sample['ewg_count']    = [r[0] for r in results]
cep_sample['ewg_weighted'] = [r[1] for r in results]
cep_sample['mol_weight']   = [r[2] for r in results]
cep_sample['n_arom_rings'] = [r[3] for r in results]
cep_sample['conj_length']  = [r[4] for r in results]
cep_sample['logp']         = [r[5] for r in results]
cep_sample = cep_sample.dropna()
print(f"CEP after feature computation: {len(cep_sample)} molecules")

# ── 6. Merge ──────────────────────────────────────────────────────
master = pd.concat([base, cep_sample], ignore_index=True)
master = master.drop_duplicates(subset=['canonical_SMILES'])
print(f"\nMaster dataset: {len(master)} molecules")

# ── 7. Save ───────────────────────────────────────────────────────
master.to_csv("master_acceptor_dataset.csv", index=False)

print("\n=== Final Summary ===")
print(f"Total molecules: {len(master)}")
print(f"\nSource breakdown:")
print(master['source_db'].value_counts())
print(f"\nHOMO range: {master['homo_ev'].min():.3f} to {master['homo_ev'].max():.3f} eV")
print(f"LUMO range: {master['lumo_ev'].min():.3f} to {master['lumo_ev'].max():.3f} eV")
print(f"Bandgap range: {master['bandgap_ev'].min():.3f} to {master['bandgap_ev'].max():.3f} eV")
print(f"EWG count range: {master['ewg_count'].min():.0f} to {master['ewg_count'].max():.0f}")
print(f"EWG weighted range: {master['ewg_weighted'].min():.2f} to {master['ewg_weighted'].max():.2f}")
print(f"\nMissing values: {master.isnull().sum().sum()}")
print("\nSaved master_acceptor_dataset.csv")