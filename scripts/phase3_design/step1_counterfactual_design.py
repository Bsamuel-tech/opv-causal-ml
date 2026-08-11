"""
Phase 3 — Step 1: Counterfactual Molecular Design Engine
=========================================================
Uses per-molecule CATE from causal forest (not population ATE)
to compute required EWG additions for each starting molecule.
Each molecule gets its own individualized causal effect estimate.

Author: Samuel Bizimana | JUNIA ISEN
Supervisor: Dr. Kekeli N'KONOU
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import os
import sys

from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

print("=== Phase 3: Counterfactual Molecular Design Engine ===")
print("Using per-molecule CATE from causal forest\n")

# ── 1. Load molecules with pre-computed CATEs ─────────────────────
# CATEs computed in R using causal_forest.R and saved to this file
df = pd.read_csv("data/processed/experimental_with_cates.csv")
print(f"Experimental molecules with CATEs: {len(df)}")
print(f"Mean CATE: {df['CATE'].mean():.6f} eV")
print(f"CATE range: {df['CATE'].min():.6f} to {df['CATE'].max():.6f} eV")

# ── 2. Filter to synthesizable starting molecules ─────────────────
def get_sa(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return 99
        return sascorer.calculateScore(mol)
    except:
        return 99

df['sascore'] = df['canonical_SMILES'].apply(get_sa)
good_starts   = df[df['sascore'] < 4.0].copy()
good_starts   = good_starts.sort_values('sascore')
print(f"\nStarting molecules with SAScore < 4.0: {len(good_starts)}")

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

# ── 4. Compute EWG features for modified molecules ────────────────
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

# ── 5. Generate candidates using per-molecule CATE ────────────────
print("\nGenerating counterfactual candidates using per-molecule CATE...")
all_candidates = []

for _, row in good_starts.iterrows():
    orig_smiles = row['canonical_SMILES']
    orig_homo   = row['homo_ev']
    orig_ewg    = row['ewg_weighted']
    cate        = row['CATE']        # per-molecule causal effect
    cate_se     = row['CATE_se']     # uncertainty

    # Compute required EWG increase for -0.2 eV HOMO shift
    target_shift       = -0.2
    required_ewg_delta = target_shift / cate if cate != 0 else None
    cyano_needed       = (required_ewg_delta / 0.66
                          if required_ewg_delta is not None else None)

    new_smiles_list = add_one_cyano(orig_smiles)
    for new_smi in new_smiles_list[:3]:
        feats = get_features(new_smi)
        if feats and feats['sascore'] < 4.0:
            delta_ewg  = feats['ewg_weighted'] - orig_ewg
            # Use per-molecule CATE for prediction
            pred_shift = cate * delta_ewg
            pred_homo  = orig_homo + pred_shift

            all_candidates.append({
                'original_smiles':    orig_smiles,
                'modified_smiles':    new_smi,
                'original_homo':      orig_homo,
                'original_ewg':       orig_ewg,
                'new_ewg':            feats['ewg_weighted'],
                'delta_ewg':          round(delta_ewg, 3),
                'cate':               round(cate, 6),
                'cate_se':            round(cate_se, 6),
                'predicted_homo':     round(pred_homo, 6),
                'predicted_shift':    round(pred_shift, 6),
                'sascore':            feats['sascore'],
                'mol_weight':         feats['mol_weight'],
                'n_cyano_added':      1
            })

# ── 6. Select top 10 by SAScore ───────────────────────────────────
candidates_df = pd.DataFrame(all_candidates)
candidates_df = candidates_df.drop_duplicates(
    subset=['modified_smiles'])
candidates_df = candidates_df.sort_values('sascore').head(10)

print(f"\nTop 10 candidates (SAScore < 4.0):")
print(candidates_df[[
    'original_homo', 'cate', 'predicted_homo',
    'predicted_shift', 'sascore', 'mol_weight'
]].to_string())

candidates_df.to_csv(
    "results/tables/phase3_top10_candidates.csv",
    index=False
)
print(f"\nSaved results/tables/phase3_top10_candidates.csv")
print(f"Total candidates: {len(candidates_df)}")
print("\nNOTE: Predicted HOMO shifts use per-molecule CATE")
print("      not the population-average DML ATE")