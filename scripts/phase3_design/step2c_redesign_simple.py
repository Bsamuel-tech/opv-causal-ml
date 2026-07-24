import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
import os, sys
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

print("=== Phase 3: Counterfactual design on simple molecule ===\n")

# ── Load simple candidates ────────────────────────────────────────
simple = pd.read_csv("results/tables/phase3_simple_candidates.csv")

# Use molecule with lowest HOMO (most interesting for design)
target_mol = simple.sort_values('homo_ev').iloc[0]
print(f"Selected molecule:")
print(f"  HOMO:        {target_mol['homo_ev']} eV")
print(f"  MW:          {target_mol['mol_weight']:.1f} g/mol")
print(f"  SAScore:     {target_mol['sascore']:.3f}")
print(f"  EWG weighted:{target_mol['ewg_weighted']}")
print(f"  SMILES:      {target_mol['canonical_SMILES'][:60]}...")

# ── SMARTS reaction: add cyano to aromatic C-H ────────────────────
cyano_rxn = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C#N')

def add_cyano_stepwise(smiles, max_steps=5):
    candidates = []
    current = smiles
    for i in range(max_steps):
        mol = Chem.MolFromSmiles(current)
        if mol is None:
            break
        products = cyano_rxn.RunReactants((mol,))
        if not products:
            break
        new_mols = []
        for prod_tuple in products:
            try:
                prod = prod_tuple[0]
                Chem.SanitizeMol(prod)
                smi = Chem.MolToSmiles(prod)
                if smi and smi != current:
                    new_mols.append(smi)
            except:
                continue
        if not new_mols:
            break
        current = new_mols[0]

        mol_new = Chem.MolFromSmiles(current)
        cyano_count = len(mol_new.GetSubstructMatches(
            Chem.MolFromSmarts('[C]#[N]')))
        carbonyl = len(mol_new.GetSubstructMatches(
            Chem.MolFromSmarts('[C]=[O]')))
        halogen = len(mol_new.GetSubstructMatches(
            Chem.MolFromSmarts('[F,Cl,Br,I]')))
        ester = len(mol_new.GetSubstructMatches(
            Chem.MolFromSmarts('[OX2][CX3](=[OX1])')))
        ewg_w = (cyano_count*0.66 + carbonyl*0.50 +
                 halogen*0.23 + ester*0.45)
        sa = sascorer.calculateScore(mol_new)
        mw = Descriptors.MolWt(mol_new)

        candidates.append({
            'n_cyano_added': i + 1,
            'smiles':        current,
            'ewg_weighted':  ewg_w,
            'mol_weight':    mw,
            'sascore':       sa,
            'sa_pass':       sa < 4.0
        })
        print(f"  +{i+1} cyano: MW={mw:.1f} SA={sa:.2f} "
              f"EWG={ewg_w:.2f} {'✓' if sa < 4.0 else '✗'}")

    return candidates

print("\nGenerating candidates...")
candidates = add_cyano_stepwise(target_mol['canonical_SMILES'])

df = pd.DataFrame(candidates)
df['original_smiles'] = target_mol['canonical_SMILES']
df['original_homo']   = target_mol['homo_ev']

passing = df[df['sa_pass']].copy()
print(f"\n{len(passing)} candidates pass SAScore < 4.0")
print(f"{len(df)} total candidates generated")

df.to_csv("results/tables/phase3_simple_all.csv", index=False)
if len(passing) > 0:
    passing.to_csv(
        "results/tables/phase3_simple_passing.csv",
        index=False
    )
    print("\nPassing candidates:")
    print(passing[['n_cyano_added','mol_weight',
                   'sascore','ewg_weighted']].to_string())
print("\nSaved phase3_simple_all.csv")