"""
Phase 3 — Step 2: xTB Quantum Chemistry Validation (FIXED)
===========================================================
Correct comparison: delta_HOMO_ML vs delta_HOMO_xTB

For each candidate:
  - Run xTB on the PARENT molecule  -> xtb_parent_homo
  - Run xTB on the MODIFIED molecule -> xtb_modified_homo
  - delta_HOMO_xTB = xtb_modified_homo - xtb_parent_homo
  - delta_HOMO_ML  = predicted_shift (CATE * delta_EWG)
  - error = |delta_HOMO_ML - delta_HOMO_xTB|

This is the correct shift-level validation. It isolates the effect
of the structural modification from the experimental-to-xTB energy
scale difference, so no calibration offset is needed for the shift
comparison.
"""

import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT           = Path(__file__).resolve().parent.parent.parent
RESULTS_TABLES = ROOT / "results" / "tables"
RESULTS_TABLES.mkdir(parents=True, exist_ok=True)

print("=== Phase 3: xTB Shift Validation (parent + modified) ===\n")

XTB_OFFSET    = 4.6197  # kept for absolute HOMO reporting only
MAE_THRESHOLD = 0.25    # eV — applied to shift error

def run_xtb_homo(smiles, mol_id, seed=42):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        result = AllChem.EmbedMolecule(mol, params)
        if result != 0:
            params2 = AllChem.ETKDG()
            params2.randomSeed = seed
            AllChem.EmbedMolecule(mol, params2)
        AllChem.MMFFOptimizeMolecule(mol)
        tmpdir   = tempfile.mkdtemp()
        xyz_file = os.path.join(tmpdir, f"mol_{mol_id}.xyz")
        conf     = mol.GetConformer()
        with open(xyz_file, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\nmolecule\n")
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")
        cmd = ["xtb", xyz_file, "--gfn", "2", "--sp", "--silent"]
        res = subprocess.run(cmd, capture_output=True, text=True,
                             encoding='utf-8', errors='ignore',
                             cwd=tmpdir, timeout=300)
        output = res.stdout + res.stderr
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
        print(f"  xTB error: {e}")
        return None

df = pd.read_csv(RESULTS_TABLES / "phase3_top10_candidates.csv")
print(f"Candidates: {len(df)}")
print(f"Offset for absolute HOMO reporting: {XTB_OFFSET} eV")
print(f"Shift MAE threshold: {MAE_THRESHOLD} eV\n")

results = []
for i, row in df.iterrows():
    n          = len(results) + 1
    orig_smi   = row['original_smiles']
    mod_smi    = row['modified_smiles']
    orig_homo  = row['original_homo']
    pred_shift = row['predicted_shift']
    cate       = row['cate']
    cate_se    = row['cate_se']

    print(f"Candidate {n}/10 (SAScore={row['sascore']:.3f})")
    print(f"  ML predicted shift: {pred_shift:+.6f} eV")
    print(f"  Running xTB on parent molecule...")
    xtb_parent_raw = run_xtb_homo(orig_smi, f"{n}_parent")

    print(f"  Running xTB on modified molecule...")
    xtb_mod_raw    = run_xtb_homo(mod_smi,  f"{n}_mod")

    if xtb_parent_raw is not None and xtb_mod_raw is not None:
        xtb_parent_cal = xtb_parent_raw + XTB_OFFSET
        xtb_mod_cal    = xtb_mod_raw    + XTB_OFFSET
        xtb_shift      = xtb_mod_cal - xtb_parent_cal   # delta_HOMO_xTB
        shift_error    = abs(pred_shift - xtb_shift)
        passed         = shift_error < MAE_THRESHOLD
        print(f"  xTB parent (cal):   {xtb_parent_cal:.4f} eV")
        print(f"  xTB modified (cal): {xtb_mod_cal:.4f} eV")
        print(f"  xTB shift:          {xtb_shift:+.4f} eV")
        print(f"  Shift error:        {shift_error:.4f} eV  {'PASS' if passed else 'FAIL'}")
    else:
        xtb_parent_cal = xtb_mod_cal = xtb_shift = shift_error = None
        passed = None
        print(f"  xTB: failed (parent={'ok' if xtb_parent_raw else 'fail'}, "
              f"mod={'ok' if xtb_mod_raw else 'fail'})")

    results.append({
        'candidate':              n,
        'original_smiles':        orig_smi,
        'modified_smiles':        mod_smi,
        'original_homo_exp_ev':   orig_homo,
        'per_molecule_cate':      cate,
        'cate_se':                cate_se,
        'ml_predicted_shift_ev':  pred_shift,
        'xtb_parent_homo_ev':     xtb_parent_cal,
        'xtb_modified_homo_ev':   xtb_mod_cal,
        'xtb_shift_ev':           xtb_shift,
        'shift_error_ev':         shift_error,
        'pass_025_shift':         passed,
        'sascore':                row['sascore'],
    })
    print()

results_df = pd.DataFrame(results)
valid      = results_df.dropna(subset=['xtb_shift_ev'])

print("=== Validation Summary ===")
print(results_df[[
    'candidate', 'ml_predicted_shift_ev',
    'xtb_shift_ev', 'shift_error_ev', 'pass_025_shift'
]].to_string())

if len(valid) > 0:
    mae    = valid['shift_error_ev'].mean()
    n_pass = int(valid['pass_025_shift'].sum())
    corr   = valid['ml_predicted_shift_ev'].corr(valid['xtb_shift_ev'])
    print(f"\nShift MAE (ML vs xTB): {mae:.4f} eV")
    print(f"Shift correlation:     {corr:.4f}")
    print(f"Candidates passing:    {n_pass}/{len(valid)}")
    print(f"Threshold:             {MAE_THRESHOLD} eV")

    # --- Zero-shift baseline comparison ---
    # Baseline: predict zero shift for every candidate (no-change assumption)
    valid = valid.copy()
    valid['baseline_error_ev'] = valid['xtb_shift_ev'].abs()
    valid['baseline_pass']     = valid['baseline_error_ev'] < MAE_THRESHOLD
    baseline_mae    = valid['baseline_error_ev'].mean()
    baseline_n_pass = int(valid['baseline_pass'].sum())
    print(f"\n--- Zero-shift baseline (predict 0 for all candidates) ---")
    print(f"Baseline MAE:          {baseline_mae:.4f} eV")
    print(f"Baseline passing:      {baseline_n_pass}/{len(valid)}")
    print(f"ML vs baseline MAE:    {mae:.4f} vs {baseline_mae:.4f} eV")
    if abs(mae - baseline_mae) < 0.01:
        print("NOTE: ML model does not outperform zero-shift baseline on this validation set.")
        print("      Results should be interpreted as exploratory.")

out_path = RESULTS_TABLES / "phase3_xtb_shift_validation.csv"
results_df.to_csv(out_path, index=False)
print(f"\nSaved {out_path}")
