import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem

print("=== xTB Validation: Top 10 Candidates ===\n")

df = pd.read_csv("results/tables/phase3_top10_candidates.csv")
print(f"Candidates to validate: {len(df)}")

# Calibration offset from Step 4b
xtb_offset = 4.6197

# CATE for each original molecule
# We use population DML ATE as approximation
dml_ate = -0.01201265

def ml_predict_homo(orig_homo, orig_ewg, new_ewg):
    delta_ewg = new_ewg - orig_ewg
    return orig_homo + (dml_ate * delta_ewg)

def run_xtb_homo(smiles, mol_id):
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
            f.write(f"{mol.GetNumAtoms()}\nmolecule\n")
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} "
                        f"{pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

        cmd = ["xtb", xyz_file, "--gfn", "2", "--sp", "--silent"]
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            encoding='utf-8', errors='ignore',
            cwd=tmpdir, timeout=180
        )
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

results = []
for i, row in df.iterrows():
    orig_homo  = row['original_homo']
    orig_ewg   = row['original_ewg']
    new_ewg    = row['new_ewg']
    sascore    = row['sascore']
    ml_homo    = ml_predict_homo(orig_homo, orig_ewg, new_ewg)
    ml_shift   = ml_homo - orig_homo

    print(f"\nCandidate {i+1}/10 (SA={sascore:.3f}):")
    print(f"  Original HOMO: {orig_homo:.3f} eV")
    print(f"  ML predicted:  {ml_homo:.4f} eV (shift={ml_shift:.4f})")
    print(f"  Running xTB...")

    xtb_raw = run_xtb_homo(row['modified_smiles'], i)

    if xtb_raw is not None:
        xtb_cal   = xtb_raw + xtb_offset
        xtb_shift = xtb_cal - orig_homo
        error     = abs(ml_homo - xtb_cal)
        passed    = error < 0.25
        print(f"  xTB calibrated: {xtb_cal:.4f} eV (shift={xtb_shift:.4f})")
        print(f"  Error: {error:.4f} eV {'✅' if passed else '❌'}")
    else:
        xtb_cal = xtb_shift = error = None
        passed  = None
        print(f"  xTB: failed")

    results.append({
        'candidate':          i + 1,
        'original_smiles':    row['original_smiles'],
        'modified_smiles':    row['modified_smiles'],
        'original_homo_ev':   orig_homo,
        'ml_predicted_homo':  ml_homo,
        'ml_homo_shift_ev':   ml_shift,
        'xtb_homo_ev':        xtb_cal,
        'xtb_homo_shift_ev':  xtb_shift,
        'sascore':            sascore,
        'abs_error_ev':       error,
        'pass_025':           passed,
        'synthesizable':      'YES' if sascore < 4.0 else 'NO'
    })

results_df = pd.DataFrame(results)
valid = results_df.dropna(subset=['xtb_homo_ev'])

print("\n=== Validation Summary ===")
print(results_df[['candidate', 'original_homo_ev', 'ml_predicted_homo',
                   'xtb_homo_ev', 'abs_error_ev',
                   'sascore', 'pass_025']].to_string())

if len(valid) > 0:
    mae = valid['abs_error_ev'].mean()
    n_pass = valid['pass_025'].sum()
    print(f"\nMAE (ML vs xTB): {mae:.4f} eV")
    print(f"Candidates passing 0.25 eV: {n_pass}/{len(valid)}")
    print(f"Direction confirmed: {(valid['xtb_homo_shift_ev'] < 0).sum()}/{len(valid)}")

results_df.to_csv("results/tables/phase3_top10_validation.csv", index=False)
print("\nSaved phase3_top10_validation.csv")