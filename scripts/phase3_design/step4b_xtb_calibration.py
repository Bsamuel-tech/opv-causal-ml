import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem

print("=== Step 4b: xTB Calibration against experimental HOMO ===\n")

# Load your experimental molecules
df = pd.read_csv("results/tables/phase3_simple_candidates.csv")
print(f"Experimental molecules for calibration: {len(df)}")

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

# Run xTB on experimental molecules
print("Running xTB on experimental molecules...")
xtb_homos = []
for i, row in df.iterrows():
    xtb = run_xtb_homo(row['canonical_SMILES'], i)
    xtb_homos.append(xtb)
    if xtb:
        print(f"  Mol {i+1}: exp={row['homo_ev']:.3f} xTB={xtb:.3f} "
              f"offset={row['homo_ev']-xtb:.3f}")

df['xtb_homo'] = xtb_homos
df = df.dropna(subset=['xtb_homo'])

if len(df) > 0:
    # Compute systematic offset
    offsets = df['homo_ev'] - df['xtb_homo']
    mean_offset = offsets.mean()
    std_offset  = offsets.std()
    print(f"\n=== Calibration Results ===")
    print(f"Mean offset (exp - xTB): {mean_offset:.4f} eV")
    print(f"Std of offset:           {std_offset:.4f} eV")
    print(f"This offset will be applied to calibrate xTB predictions")

    # Apply calibration to Phase 3 candidates
    candidates = pd.read_csv("results/tables/phase3_xtb_validation.csv")
    candidates['xtb_homo_calibrated'] = candidates['xtb_homo'] + mean_offset
    candidates['calibrated_error'] = abs(
        candidates['ml_predicted_homo'] - candidates['xtb_homo_calibrated']
    )
    candidates['pass_calibrated'] = candidates['calibrated_error'] < 0.25

    print(f"\n=== Calibrated Validation Results ===")
    print(candidates[['n_cyano_added', 'ml_predicted_homo',
                       'xtb_homo', 'xtb_homo_calibrated',
                       'calibrated_error', 'pass_calibrated']].to_string())

    mae_cal = candidates['calibrated_error'].mean()
    print(f"\nCalibrated MAE: {mae_cal:.4f} eV")
    print(f"Target MAE:     0.25 eV")
    print(f"Validation:     {'PASSED' if mae_cal < 0.25 else 'FAILED'}")

    # Save
    candidates.to_csv("results/tables/phase3_xtb_calibrated.csv", index=False)
    print("\nSaved phase3_xtb_calibrated.csv")