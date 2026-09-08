import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski
from sklearn.feature_selection import SelectKBest, f_regression

df = pd.read_csv("data/processed/master_acceptor_dataset.csv")
df = df[df['ewg_count'] > 0].copy().reset_index(drop=True)
print("Molecules:", len(df))

tpsa_list    = []
hbd_list     = []
hba_list     = []
rot_list     = []
heavy_list   = []
hetero_list  = []
rings_list   = []
molar_list   = []
bertz_list   = []
valid_idx    = []

for i, smi in enumerate(df['canonical_SMILES']):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            continue
        tpsa_list.append(Descriptors.TPSA(mol))
        hbd_list.append(Lipinski.NumHDonors(mol))
        hba_list.append(Lipinski.NumHAcceptors(mol))
        rot_list.append(Descriptors.NumRotatableBonds(mol))
        heavy_list.append(Descriptors.HeavyAtomCount(mol))
        hetero_list.append(Descriptors.NumHeteroatoms(mol))
        rings_list.append(rdMolDescriptors.CalcNumRings(mol))
        molar_list.append(Descriptors.MolMR(mol))
        bertz_list.append(Descriptors.BertzCT(mol))
        valid_idx.append(i)
    except Exception as e:
        print(f"Row {i} error: {e}")
        continue

print("Valid molecules:", len(valid_idx))

feat_df = pd.DataFrame({
    'tpsa':            tpsa_list,
    'hbd':             hbd_list,
    'hba':             hba_list,
    'rotatable_bonds': rot_list,
    'heavy_atoms':     heavy_list,
    'heteroatoms':     hetero_list,
    'fused_rings':     rings_list,
    'molar_ref':       molar_list,
    'bertz_ct':        bertz_list,
})

df_valid = df.iloc[valid_idx].reset_index(drop=True)

print("Features computed:", feat_df.shape[1])
print(feat_df.describe().round(2))

X = feat_df.values
y = df_valid['homo_ev'].values

selector = SelectKBest(score_func=f_regression, k=5)
selector.fit(X, y)
selected = feat_df.columns[selector.get_support()].tolist()
scores   = selector.scores_[selector.get_support()]

print("\n=== Top 5 features most informative for HOMO ===")
for feat, score in sorted(zip(selected, scores),
                          key=lambda x: x[1], reverse=True):
    print(f"  {feat:<20} F-score: {score:.2f}")

feat_df['canonical_SMILES'] = df_valid['canonical_SMILES'].values
feat_df['homo_ev']          = df_valid['homo_ev'].values
feat_df['ewg_weighted']     = df_valid['ewg_weighted'].values
feat_df.to_csv("results/tables/extended_features.csv", index=False)
print("\nSaved results/tables/extended_features.csv")
print("\nRecommended additional covariates for causal forest:")
for f in selected:
    print(f"  {f}")