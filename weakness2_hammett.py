import pandas as pd
import numpy as np
from rdkit import Chem
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler

df = pd.read_csv("data/processed/master_acceptor_dataset.csv")
df = df[df['ewg_count'] > 0].copy()
print("Molecules:", len(df))

# Count each EWG type per molecule
def count_ewgs(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return [0, 0, 0, 0, 0, 0]
    cyano    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]#[N]')))
    carbonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[C]=[O]')))
    nitro    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[N+](=O)[O-]')))
    halogen  = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F,Cl,Br,I]')))
    ester    = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2][CX3](=[OX1])')))
    sulfonyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[SX4](=O)(=O)')))
    return [cyano, carbonyl, nitro, halogen, ester, sulfonyl]

print("Computing EWG counts...")
ewg_cols = ['n_cyano','n_carbonyl','n_nitro','n_halogen','n_ester','n_sulfonyl']
ewg_matrix = df['canonical_SMILES'].apply(count_ewgs)
df[ewg_cols] = pd.DataFrame(ewg_matrix.tolist(), index=df.index)

print("EWG counts computed:")
print(df[ewg_cols].describe().round(3))

# LASSO regression: predict LUMO from EWG counts + structural controls
controls = ['mol_weight', 'n_arom_rings', 'conj_length']
features = ewg_cols + controls

X = df[features].values
y = df['lumo_ev'].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

lasso = LassoCV(cv=5, random_state=42, max_iter=5000)
lasso.fit(X_scaled, y)

# Extract empirical EWG weights (unscaled for EWG columns only)
empirical = {}
for i, col in enumerate(ewg_cols):
    empirical[col] = round(lasso.coef_[i], 4)

hammett = {
    'n_cyano':    0.66,
    'n_carbonyl': 0.50,
    'n_nitro':    0.78,
    'n_halogen':  0.23,
    'n_ester':    0.45,
    'n_sulfonyl': 0.72
}

print("\n=== Hammett vs Empirical EWG weights ===")
print(f"{'EWG':<12} {'Hammett':>10} {'Empirical':>12}")
print("-" * 36)
emp_vals = []
ham_vals = []
for key in ewg_cols:
    h = hammett[key]
    e = empirical[key]
    emp_vals.append(e)
    ham_vals.append(h)
    print(f"{key:<12} {h:>10.3f} {e:>12.4f}")

# Correlation between Hammett and empirical
corr = np.corrcoef(ham_vals, emp_vals)[0, 1]
print(f"\nCorrelation (Hammett vs Empirical): {corr:.4f}")

if corr > 0.8:
    print("Result: Hammett values VALIDATED by data (correlation > 0.8)")
    print("Conclusion: Keep Hammett weights in manuscript")
elif corr > 0.5:
    print("Result: Moderate correlation - Hammett values partially supported")
    print("Conclusion: Report both and discuss discrepancy")
else:
    print("Result: Low correlation - Hammett values not supported by this dataset")
    print("Conclusion: Consider using empirical weights")

# Save results
results = pd.DataFrame({
    'ewg_group':       ewg_cols,
    'hammett_sigma':   [hammett[k] for k in ewg_cols],
    'empirical_coef':  [empirical[k] for k in ewg_cols],
    'correlation':     [corr] * 6
})
results.to_csv("results/tables/hammett_validation.csv", index=False)
print("\nSaved results/tables/hammett_validation.csv")
print(results)