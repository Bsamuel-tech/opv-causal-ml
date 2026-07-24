import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import os, sys
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

print("=== Finding a simpler test molecule ===\n")

df = pd.read_csv("data/corrected_pipeline/master_acceptor_dataset.csv")
df = df[df['ewg_count'] > 0].copy()
df = df[df['source_db'] == 'your_data'].copy()

print(f"Experimental molecules: {len(df)}")

# Compute SAScore for all experimental molecules
def get_sa(smi):
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return 99
        return sascorer.calculateScore(mol)
    except:
        return 99

print("Computing SAScores...")
df['sascore'] = df['canonical_SMILES'].apply(get_sa)
df['mol_weight'] = df['canonical_SMILES'].apply(
    lambda s: Descriptors.MolWt(Chem.MolFromSmiles(s))
    if Chem.MolFromSmiles(s) else 99
)

# Find molecules with good SAScore and moderate size
good = df[
    (df['sascore'] < 4.0) &
    (df['mol_weight'] < 800)
].sort_values('sascore')

print(f"\nMolecules with SAScore < 4.0 and MW < 800: {len(good)}")

if len(good) > 0:
    print("\nTop 5 most synthesizable:")
    print(good[['canonical_SMILES', 'homo_ev', 'ewg_weighted',
                'mol_weight', 'sascore']].head())
    good.head(10).to_csv(
        "results/tables/phase3_simple_candidates.csv",
        index=False
    )
    print("\nSaved phase3_simple_candidates.csv")
else:
    print("Relaxing to SAScore < 5.0...")
    good = df[df['sascore'] < 5.0].sort_values('sascore')
    print(f"Found {len(good)} molecules")
    print(good[['canonical_SMILES', 'homo_ev', 'ewg_weighted',
                'mol_weight', 'sascore']].head())
    good.head(10).to_csv(
        "results/tables/phase3_simple_candidates.csv",
        index=False
    )