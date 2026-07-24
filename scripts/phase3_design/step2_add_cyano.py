import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticRings

print("=== Step 2: Adding cyano groups via SMARTS reactions ===\n")

# ── Load test molecule ────────────────────────────────────────────
info = pd.read_csv("results/tables/phase3_test_molecule.csv")
original_smiles = info['smiles'].iloc[0]
original_homo   = info['homo_ev'].iloc[0]
cyano_needed    = int(round(info['cyano_needed'].iloc[0]))

print(f"Original SMILES: {original_smiles[:60]}...")
print(f"Original HOMO:   {original_homo} eV")
print(f"Cyano groups to add: {cyano_needed}")

# ── SMARTS reaction: add cyano to aromatic C-H ────────────────────
# This replaces an aromatic C-H bond with C-C#N
cyano_rxn = AllChem.ReactionFromSmarts(
    '[cH:1]>>[c:1]C#N'
)

def add_cyano_groups(smiles, n_groups):
    """Add up to n_groups cyano groups to available aromatic C-H positions."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    generated = []
    current_smiles = smiles

    for i in range(min(n_groups, 5)):  # max 5 additions
        mol = Chem.MolFromSmiles(current_smiles)
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
                if smi and smi != current_smiles:
                    new_mols.append(smi)
            except:
                continue

        if not new_mols:
            break

        # Take first valid product and continue
        current_smiles = new_mols[0]
        generated.append({
            'step':         i + 1,
            'smiles':       current_smiles,
            'n_cyano_added': i + 1
        })
        print(f"  Step {i+1}: Added cyano group {i+1}")

    return generated

print("\nGenerating counterfactual molecules...")
candidates = add_cyano_groups(original_smiles, cyano_needed)

if candidates:
    print(f"\nGenerated {len(candidates)} candidate molecules")
    df = pd.DataFrame(candidates)

    # Compute basic properties for each
    def get_props(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None, None
        return Descriptors.MolWt(mol), CalcNumAromaticRings(mol)

    df[['mol_weight', 'n_arom_rings']] = df['smiles'].apply(
        lambda s: pd.Series(get_props(s))
    )
    df['original_homo'] = original_homo
    df.to_csv("results/tables/phase3_candidates.csv", index=False)
    print("\nCandidate molecules:")
    print(df[['step', 'n_cyano_added', 'mol_weight']].to_string())
    print("\nSaved results/tables/phase3_candidates.csv")
else:
    print("No candidates generated — molecule may have no free aromatic C-H positions")
    print("This is common for highly substituted acceptor molecules.")
    print("We will try a different molecule from the dataset.")