import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticRings

print("=== Step 3: SAScore synthesizability filtering ===\n")

# ── Load SA Score function from RDKit ─────────────────────────────
from rdkit.Chem import RDConfig
import os, sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

# ── Load candidates ───────────────────────────────────────────────
df = pd.read_csv("results/tables/phase3_candidates.csv")
original_homo = df['original_homo'].iloc[0]
print(f"Candidates to evaluate: {len(df)}")

# ── Compute SAScore and EWG features for each candidate ───────────
def compute_features(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None

        sa = sascorer.calculateScore(mol)

        # EWG features
        cyano    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]#[N]')))
        carbonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=[O]')))
        nitro    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[N+](=O)[O-]')))
        halogen  = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F,Cl,Br,I]')))
        ester    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2][CX3](=[OX1])')))
        sulfonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[SX4](=[OX1])(=[OX1])')))

        ewg_weighted = (cyano*0.66 + carbonyl*0.50 + nitro*0.78 +
                       halogen*0.23 + ester*0.45 + sulfonyl*0.72)

        mol_weight   = Descriptors.MolWt(mol)
        n_arom_rings = CalcNumAromaticRings(mol)
        conj_length  = rdMolDescriptors.CalcNumRotatableBonds(mol)
        logp         = Descriptors.MolLogP(mol)

        return {
            'sascore':      sa,
            'ewg_weighted': ewg_weighted,
            'mol_weight':   mol_weight,
            'n_arom_rings': n_arom_rings,
            'conj_length':  conj_length,
            'logp':         logp
        }
    except:
        return None

print("Computing SAScore and molecular features...")
results = []
for _, row in df.iterrows():
    feats = compute_features(row['smiles'])
    if feats:
        results.append({
            'step':         row['step'],
            'n_cyano_added': row['n_cyano_added'],
            'smiles':       row['smiles'],
            'original_homo': original_homo,
            **feats
        })

results_df = pd.DataFrame(results)

print("\n=== SAScore Results ===")
print("(SAScore < 4.0 = synthesizable, < 6.0 = moderate difficulty)")
for _, r in results_df.iterrows():
    flag = "✓ PASS" if r['sascore'] < 6.0 else "✗ FAIL"
    print(f"  +{int(r['n_cyano_added'])} cyano: SAScore={r['sascore']:.2f} "
          f"EWG={r['ewg_weighted']:.2f} {flag}")

# Filter
passing = results_df[results_df['sascore'] < 6.0].copy()
print(f"\n{len(passing)} molecules pass SAScore < 6.0 threshold")

# Save all and passing
results_df.to_csv("results/tables/phase3_all_candidates.csv", index=False)
passing.to_csv("results/tables/phase3_passing_candidates.csv", index=False)
print("Saved phase3_all_candidates.csv and phase3_passing_candidates.csv")