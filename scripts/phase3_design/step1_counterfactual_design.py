"""
Phase 3 — Step 1: Counterfactual Molecular Design Engine
=========================================================
Given synthesizable organic acceptor molecules, this script:
1. Finds molecules with SAScore < 4.0 from experimental dataset
2. Adds cyano groups via SMARTS reaction [cH]>>[c]C#N
3. Filters candidates by SAScore < 4.0
4. Saves top 10 candidates for xTB validation

Note: SYBA (Synthesizability by Bayesian Approach) is not used
because the package is no longer publicly available as of August 2026.
SAScore is used as the sole synthesizability filter, consistent with
the majority of recent molecular design literature.

Author: Samuel Bizimana | JUNIA ISEN
Supervisor: Dr. Kekeli N'KONOU
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import os
import sys

# SAScore
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

print("=== Phase 3: Counterfactual Molecular Design Engine ===\n")

# ── 1. Load experimental molecules ───────────────────────────────
df = pd.read_csv("data/corrected_pipeline/master_acceptor_dataset.csv")
df = df[df['source_db'] == 'your_data'].copy()
print(f"Experimental molecules: {len(df)}")

# ── 2. Compute SAScore for starting molecules ─────────────────────
def get_sa(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return 99
        return sascorer.calculateScore(mol)
    except:
        return 99

df['sascore'] = df['canonical_SMILES'].apply(get_sa)
good_starts   = df[df['sascore'] < 4.0].sort_values('sascore')
print(f"Starting molecules with SAScore < 4.0: {len(good_starts)}")

# ── 3. SMARTS reaction: add cyano to aromatic C-H ─────────────────
cyano_rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C#N')

def add_one_cyano(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    products = cyano_rxn.RunReactants((mol,))
    results  = []
    for prod_tuple in products:
        try:
            prod = prod_tuple[0]
            Chem.SanitizeMol(prod)
            smi = Chem.MolToSmiles(prod)
            if smi and smi != smiles:
                results.append(smi)
        except:
            continue
    return list(set(results))

# ── 4. Compute features for each candidate ────────────────────────
def get_features(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        cyano    = len(mol.GetSubstructMatches(
                       Chem.MolFromSmarts('[C]#[N]')))
        carbonyl = len(mol.GetSubstructMatches(
                       Chem.MolFromSmarts('[C]=[O]')))
        halogen  = len(mol.GetSubstructMatches(
                       Chem.MolFromSmarts('[F,Cl,Br,I]')))
        ester    = len(mol.GetSubstructMatches(
                       Chem.MolFromSmarts('[OX2][CX3](=[OX1])')))
        ewg_w    = (cyano*0.66 + carbonyl*0.50 +
                    halogen*0.23 + ester*0.45)
        sa       = sascorer.calculateScore(mol)
        mw       = Descriptors.MolWt(mol)
        return {
            'ewg_weighted': round(ewg_w, 3),
            'sascore':      round(sa, 3),
            'mol_weight':   round(mw, 3)
        }
    except:
        return None

# ── 5. Generate candidates ────────────────────────────────────────
print("\nGenerating counterfactual candidates...")
all_candidates = []

for _, row in good_starts.iterrows():
    orig_smiles = row['canonical_SMILES']
    orig_homo   = row['homo_ev']
    orig_ewg    = row['ewg_weighted']

    new_smiles_list = add_one_cyano(orig_smiles)
    for new_smi in new_smiles_list[:3]:
        feats = get_features(new_smi)
        if feats and feats['sascore'] < 4.0:
            all_candidates.append({
                'original_smiles': orig_smiles,
                'modified_smiles': new_smi,
                'original_homo':   orig_homo,
                'original_ewg':    orig_ewg,
                'new_ewg':         feats['ewg_weighted'],
                'sascore':         feats['sascore'],
                'mol_weight':      feats['mol_weight'],
                'n_cyano_added':   1
            })

# ── 6. Select top 10 by SAScore ───────────────────────────────────
candidates_df = pd.DataFrame(all_candidates)
candidates_df = candidates_df.drop_duplicates(
    subset=['modified_smiles'])
candidates_df = candidates_df.sort_values('sascore').head(10)

print(f"\nTop 10 candidates (SAScore < 4.0):")
print(candidates_df[[
    'original_homo', 'new_ewg',
    'sascore', 'mol_weight'
]].to_string())

candidates_df.to_csv(
    "results/tables/phase3_top10_candidates.csv",
    index=False
)
print(f"\nSaved results/tables/phase3_top10_candidates.csv")
print(f"Total candidates: {len(candidates_df)}")