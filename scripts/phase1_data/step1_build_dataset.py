import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
from rdkit.Chem.Scaffolds import MurckoScaffold

print("=== Step 1: Building correct acceptor dataset ===\n")

# ── 1. Load your 1,573-row dataset ────────────────────────────────
df = pd.read_excel(
    r"C:\Users\Samuel Bizimana\OneDrive\Desktop\Research Training\data\Data_Merged_with_SMILES.xlsx"
)
print(f"Loaded: {len(df)} rows")

# ── 2. Keep acceptor columns only ─────────────────────────────────
acc = df[['SMILES_acc', 'HOMO_A', 'LUMO_A', 'EgA_opt']].copy()
acc = acc.rename(columns={
    'SMILES_acc': 'SMILES',
    'HOMO_A':     'homo_ev',
    'LUMO_A':     'lumo_ev',
    'EgA_opt':    'bandgap_ev'
})
acc['source_db'] = 'your_data'
print(f"After selecting acceptor columns: {len(acc)} rows")

# ── 3. Validate and canonicalize SMILES ───────────────────────────
def canonical(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        return Chem.MolToSmiles(mol) if mol else None
    except:
        return None

acc['canonical_SMILES'] = acc['SMILES'].apply(canonical)
before = len(acc)
acc = acc.dropna(subset=['canonical_SMILES'])
print(f"After SMILES validation: {len(acc)} rows (dropped {before - len(acc)})")

# ── 4. Deduplicate by canonical SMILES ────────────────────────────
acc = acc.drop_duplicates(subset=['canonical_SMILES'])
print(f"After deduplication: {len(acc)} unique acceptor molecules")

# ── 5. Physical range filter ──────────────────────────────────────
acc = acc[(acc['homo_ev'] > -8) & (acc['homo_ev'] < -3)]
acc = acc[(acc['lumo_ev'] > -6) & (acc['lumo_ev'] < -1)]
acc = acc[(acc['bandgap_ev'] > 0.5) & (acc['bandgap_ev'] < 4.0)]
print(f"After physical range filter: {len(acc)} rows")

# ── 6. Compute EWG features (your treatment variables) ────────────
def compute_ewg_features(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None, None, None, None, None, None

        # EWG count: cyano, carbonyl, nitro, halogen, ester, sulfonyl
        cyano    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]#[N]')))
        carbonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=[O]')))
        nitro    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[N+](=O)[O-]')))
        halogen  = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F,Cl,Br,I]')))
        ester    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2][CX3](=[OX1])')))
        sulfonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[SX4](=[OX1])(=[OX1])')))
        ewg_count = cyano + carbonyl + nitro + halogen + ester + sulfonyl

        # Hammett-weighted EWG score (sigma_para values from literature)
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

        # Structural confounders
        mol_weight   = Descriptors.MolWt(mol)
        n_arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        conj_length  = rdMolDescriptors.CalcNumRotatableBonds(mol)
        logp         = Descriptors.MolLogP(mol)

        return ewg_count, ewg_weighted, mol_weight, n_arom_rings, conj_length, logp
    except:
        return None, None, None, None, None, None

print("\nComputing EWG features for your molecules...")
results = acc['canonical_SMILES'].apply(compute_ewg_features)
acc['ewg_count']    = [r[0] for r in results]
acc['ewg_weighted'] = [r[1] for r in results]
acc['mol_weight']   = [r[2] for r in results]
acc['n_arom_rings'] = [r[3] for r in results]
acc['conj_length']  = [r[4] for r in results]
acc['logp']         = [r[5] for r in results]
acc = acc.dropna()
print(f"After feature computation: {len(acc)} rows")

# ── 7. Summary ─────────────────────────────────────────────────────
print("\n=== Dataset Summary ===")
print(f"Total molecules: {len(acc)}")
print(f"\nHOMO range: {acc['homo_ev'].min():.3f} to {acc['homo_ev'].max():.3f} eV")
print(f"LUMO range: {acc['lumo_ev'].min():.3f} to {acc['lumo_ev'].max():.3f} eV")
print(f"Bandgap range: {acc['bandgap_ev'].min():.3f} to {acc['bandgap_ev'].max():.3f} eV")
print(f"\nEWG count range: {acc['ewg_count'].min():.0f} to {acc['ewg_count'].max():.0f}")
print(f"EWG weighted range: {acc['ewg_weighted'].min():.2f} to {acc['ewg_weighted'].max():.2f}")
print(f"\nSource breakdown:")
print(acc['source_db'].value_counts())

# ── 8. Save ───────────────────────────────────────────────────────
acc.to_csv("acceptor_base.csv", index=False)
print(f"\nSaved acceptor_base.csv with {len(acc)} molecules")
print("Columns:", acc.columns.tolist())