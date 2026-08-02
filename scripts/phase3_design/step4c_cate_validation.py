import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem

print("=== Step 4c: xTB Validation using CATE-based ML prediction ===\n")

# The CATE for our test molecule from the causal forest was -0.027865 eV
# Original molecule: HOMO = -6.20, EWG weighted = 1.32
original_homo = -6.20
original_ewg  = 1.32
cate          = -0.027865  # molecule-specific causal effect

def ml_predict_homo_cate(ewg_weighted):
    delta_ewg = ewg_weighted - original_ewg
    return original_homo + (cate * delta_ewg)

# Calibration offset from Step 4b
xtb_offset = 4.6197

def run_xtb_homo(smiles, mol_id):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
        tmpdir = tempfile.mkdtemp()
        xyz_file = os.path.join(tmpdir, f"mol_{mol_id}.xyz")
        conf = mol.GetConformer()
        with open(xyz_file, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\nmolecule\n")
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} "
                        f"{pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
        cmd = ["xtb", xyz_file, "--gfn", "2", "--sp", "--silent"]
        result = subprocess.run(cmd, capture_output=True, text=True,
                               encoding='utf-8', errors='ignore',
                               cwd=tmpdir, timeout=120)
        output = result.stdout + result.stderr
        for line in output.split('\n'):
            if '(HOMO)' in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == '(HOMO)':
                        try:
                            return float(parts[i - 1])
                        except:
                            continue
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None

# Load candidates
df = pd.read_csv("results/tables/phase3_simple_passing.csv")
print(f"Candidates: {len(df)}")
print(f"Using CATE = {cate} eV per unit EWG score")
print(f"Using xTB calibration offset = {xtb_offset} eV\n")

results = []
for i, row in df.iterrows():
    n_cyano  = int(row['n_cyano_added'])
    smiles   = row['smiles']
    ewg_w    = row['ewg_weighted']
    ml_homo  = ml_predict_homo_cate(ewg_w)

    print(f"Candidate +{n_cyano} cyano:")
    print(f"  EWG weighted:       {ewg_w:.2f}")
    print(f"  CATE-based ML HOMO: {ml_homo:.4f} eV")

    xtb_raw = run_xtb_homo(smiles, n_cyano)
    if xtb_raw is not None:
        xtb_cal = xtb_raw + xtb_offset
        error   = abs(ml_homo - xtb_cal)
        print(f"  xTB raw HOMO:       {xtb_raw:.4f} eV")
        print(f"  xTB calibrated:     {xtb_cal:.4f} eV")
        print(f"  Absolute error:     {error:.4f} eV "
              f"{'✅' if error < 0.25 else '❌'}")
    else:
        xtb_cal = None
        error   = None
        print(f"  xTB: failed")

    results.append({
        'n_cyano_added':    n_cyano,
        'ewg_weighted':     ewg_w,
        'ml_homo_cate':     ml_homo,
        'xtb_homo_raw':     xtb_raw,
        'xtb_homo_cal':     xtb_cal,
        'abs_error':        error,
        'pass_025':         error < 0.25 if error is not None else None
    })

results_df = pd.DataFrame(results)
valid = results_df.dropna(subset=['xtb_homo_cal'])

print("\n=== Final Validation Summary ===")
print(results_df[['n_cyano_added','ml_homo_cate',
                   'xtb_homo_cal','abs_error','pass_025']].to_string())

if len(valid) > 0:
    mae = valid['abs_error'].mean()
    print(f"\nMAE (CATE-ML vs calibrated xTB): {mae:.4f} eV")
    print(f"Target MAE: 0.25 eV")
    print(f"Validation: {'PASSED ✅' if mae < 0.25 else 'FAILED ❌'}")

results_df.to_csv("results/tables/phase3_xtb_cate_validation.csv", index=False)
print("\nSaved phase3_xtb_cate_validation.csv")