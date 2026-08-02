import subprocess, os, tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd

df = pd.read_csv("results/tables/phase3_simple_passing.csv")
smiles = df['smiles'].iloc[0]

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol)

tmpdir = tempfile.mkdtemp()
xyz_file = os.path.join(tmpdir, "mol.xyz")
conf = mol.GetConformer()

with open(xyz_file, 'w') as f:
    f.write(f"{mol.GetNumAtoms()}\nmolecule\n")
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

cmd = ["xtb", xyz_file, "--gfn", "2", "--sp", "--silent"]
result = subprocess.run(cmd, capture_output=True, text=True,
                        encoding='utf-8', errors='ignore',
                        cwd=tmpdir, timeout=120)

output = result.stdout + result.stderr
print("=== Lines containing HOMO or orbital energies ===")
for line in output.split('\n'):
    if any(k in line for k in ['HOMO', 'homo', 'orbital', 'eV', 'Eh']):
        print(repr(line))