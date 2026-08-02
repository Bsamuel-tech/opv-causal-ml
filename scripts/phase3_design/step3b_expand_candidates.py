import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import os, sys
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

print("=== Step 3b: Expanding to top 10 candidates ===\n")

df = pd.read_csv("results/tables/phase3_simple_candidates.csv")
print(f"Starting molecules: {len(df)}")

cyano_rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C#N')

def add_one_cyano(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    products = cyano_rxn.RunReactants((mol,))
    results = []
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
        ewg_w = (cyano*0.66 + carbonyl*0.50 +
                 halogen*0.23 + ester*0.45)
        sa    = sascorer.calculateScore(mol)
        mw    = Descriptors.MolWt(mol)
        return {'ewg_weighted': ewg_w, 'sascore': sa,
                'mol_weight': mw}
    except:
        return None

all_candidates = []

for _, row in df.iterrows():
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
                'sascore':         round(feats['sascore'], 3),
                'mol_weight':      round(feats['mol_weight'], 3),
                'n_cyano_added':   1
            })

candidates_df = pd.DataFrame(all_candidates)
candidates_df = candidates_df.drop_duplicates(subset=['modified_smiles'])
candidates_df = candidates_df.sort_values('sascore').head(10)

print(f"Top 10 candidates passing SAScore < 4.0:")
print(candidates_df[['original_homo','new_ewg',
                      'sascore','mol_weight']].to_string())

candidates_df.to_csv(
    "results/tables/phase3_top10_candidates.csv",
    index=False)
print(f"\nTotal candidates: {len(candidates_df)}")
print("Saved phase3_top10_candidates.csv")