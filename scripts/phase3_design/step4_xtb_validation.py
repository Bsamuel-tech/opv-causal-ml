import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

print("=== Step 4: xTB Validation of Phase 3 Candidates ===\n")

# ── Load passing candidates ───────────────────────────────────────
df = pd.read_csv("results/tables/phase3_simple_passing.csv")
print(f"Candidates to validate: {len(df)}")

# ── ML prediction function ────────────────────────────────────────
original_homo  = -6.20
original_ewg   = 1.32
dml_ate        = -0.01201265

def ml_predict_homo(ewg_weighted):
    delta_ewg = ewg_weighted - original_ewg
    return original_homo + (dml_ate * delta_ewg)

# ── xTB HOMO extraction ───────────────────────────────────────────
def run_xtb(smiles, mol_id):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result != 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)

        tmpdir = tempfile.mkdtemp()
        xyz_file = os.path.join(tmpdir, f"mol_{mol_id}.xyz")
        conf = mol.GetConformer()

        with open(xyz_file, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\nmolecule {mol_id}\n")
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} "
                        f"{pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

        cmd = ["xtb", xyz_file, "--gfn", "2", "--sp", "--silent"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            cwd=tmpdir,
            timeout=120
        )

        output = result.stdout + result.stderr

        # Parse HOMO in eV — line looks like:
        # 81  2.0000  -0.3948957  -10.7457 (HOMO)
        homo_ev = None
        for line in output.split('\n'):
            if '(HOMO)' in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == '(HOMO)':
                        try:
                            homo_ev = float(parts[i - 1])
                            break
                        except:
                            continue
                if homo_ev is not None:
                    break

        return homo_ev

    except Exception as e:
        print(f"  xTB error for mol {mol_id}: {e}")
        return None

# ── Run validation ────────────────────────────────────────────────
results = []

for i, row in df.iterrows():
    n_cyano  = int(row['n_cyano_added'])
    smiles   = row['smiles']
    ewg_w    = row['ewg_weighted']
    ml_homo  = ml_predict_homo(ewg_w)

    print(f"\nCandidate +{n_cyano} cyano:")
    print(f"  EWG weighted:      {ewg_w:.2f}")
    print(f"  ML predicted HOMO: {ml_homo:.4f} eV")
    print(f"  Running xTB...")

    xtb_homo = run_xtb(smiles, n_cyano)

    if xtb_homo is not None:
        error = abs(ml_homo - xtb_homo)
        print(f"  xTB HOMO:          {xtb_homo:.4f} eV")
        print(f"  Absolute error:    {error:.4f} eV")
    else:
        error = None
        print(f"  xTB: could not parse HOMO")

    results.append({
        'n_cyano_added':     n_cyano,
        'smiles':            smiles,
        'ewg_weighted':      ewg_w,
        'ml_predicted_homo': ml_homo,
        'xtb_homo':          xtb_homo,
        'abs_error':         error,
        'pass_025':          error < 0.25 if error is not None else None
    })

# ── Summary ───────────────────────────────────────────────────────
results_df = pd.DataFrame(results)
valid = results_df.dropna(subset=['xtb_homo'])

print("\n=== Validation Summary ===")
print(results_df[['n_cyano_added', 'ml_predicted_homo',
                   'xtb_homo', 'abs_error', 'pass_025']].to_string())

if len(valid) > 0:
    mae = valid['abs_error'].mean()
    print(f"\nMAE (ML vs xTB): {mae:.4f} eV")
    print(f"Target MAE:      0.25 eV")
    print(f"Validation:      {'PASSED' if mae < 0.25 else 'FAILED'}")

results_df.to_csv("results/tables/phase3_xtb_validation.csv", index=False)
print("\nSaved results/tables/phase3_xtb_validation.csv")