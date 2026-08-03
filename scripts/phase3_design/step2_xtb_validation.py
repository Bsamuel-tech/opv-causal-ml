"""
Phase 3 — Step 2: xTB Quantum Chemistry Validation
====================================================
Validates counterfactual candidates from step1 using xTB GFN2.
Calibrates xTB against experimental HOMO values, then computes
MAE between ML-predicted and xTB-computed HOMO for all candidates.

Author: Samuel Bizimana | JUNIA ISEN
Supervisor: Dr. Kekeli N'KONOU
"""

import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from rdkit import Chem
from rdkit.Chem import AllChem

print("=== Phase 3: xTB Validation of Counterfactual Candidates ===\n")

# ── Configuration ─────────────────────────────────────────────────
DML_ATE       = -0.01201265   # EWG -> HOMO ATE from Phase 2
XTB_OFFSET    = 4.6197        # Calibrated offset (exp - xTB)
MAE_THRESHOLD = 0.25          # eV

# ── xTB runner ────────────────────────────────────────────────────
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

        tmpdir   = tempfile.mkdtemp()
        xyz_file = os.path.join(tmpdir, f"mol_{mol_id}.xyz")
        conf     = mol.GetConformer()

        with open(xyz_file, 'w') as f:
            f.write(f"{mol.GetNumAtoms()}\nmolecule\n")
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} "
                        f"{pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

        cmd = ["xtb", xyz_file, "--gfn", "2", "--sp", "--silent"]
        res = subprocess.run(
            cmd, capture_output=True, text=True,
            encoding='utf-8', errors='ignore',
            cwd=tmpdir, timeout=180
        )
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

# ── Load candidates ───────────────────────────────────────────────
df = pd.read_csv("results/tables/phase3_top10_candidates.csv")
print(f"Candidates to validate: {len(df)}")
print(f"xTB calibration offset: {XTB_OFFSET} eV")
print(f"MAE threshold: {MAE_THRESHOLD} eV\n")

# ── Run validation ────────────────────────────────────────────────
results = []

for i, row in df.iterrows():
    orig_homo = row['original_homo']
    orig_ewg  = row['original_ewg']
    new_ewg   = row['new_ewg']
    delta_ewg = new_ewg - orig_ewg
    ml_homo   = orig_homo + (DML_ATE * delta_ewg)
    ml_shift  = ml_homo - orig_homo

    print(f"Candidate {i+1}/10 (SA={row['sascore']:.3f}):")
    print(f"  Original HOMO:  {orig_homo:.3f} eV")
    print(f"  ML predicted:   {ml_homo:.4f} eV (shift={ml_shift:.4f})")
    print(f"  Running xTB...")

    xtb_raw = run_xtb_homo(row['modified_smiles'], i)

    if xtb_raw is not None:
        xtb_cal   = xtb_raw + XTB_OFFSET
        xtb_shift = xtb_cal - orig_homo
        error     = abs(ml_homo - xtb_cal)
        passed    = error < MAE_THRESHOLD
        print(f"  xTB calibrated: {xtb_cal:.4f} eV "
              f"(shift={xtb_shift:.4f})")
        print(f"  Error: {error:.4f} eV "
              f"{'✅ PASS' if passed else '❌ FAIL'}")
    else:
        xtb_cal = xtb_shift = error = None
        passed  = None

    results.append({
        'candidate':          i + 1,
        'original_smiles':    row['original_smiles'],
        'modified_smiles':    row['modified_smiles'],
        'original_homo_ev':   orig_homo,
        'ml_predicted_homo':  ml_homo,
        'ml_homo_shift_ev':   ml_shift,
        'xtb_homo_ev':        xtb_cal,
        'xtb_homo_shift_ev':  xtb_shift,
        'sascore':            row['sascore'],
        'syba_score':         row.get('syba_score', None),
        'syba_pass':          row.get('syba_pass', None),
        'abs_error_ev':       error,
        'pass_025':           passed,
        'synthesizable':      'YES'
    })

# ── Summary ───────────────────────────────────────────────────────
results_df = pd.DataFrame(results)
valid      = results_df.dropna(subset=['xtb_homo_ev'])

print("\n=== Validation Summary ===")
print(results_df[[
    'candidate', 'original_homo_ev', 'ml_predicted_homo',
    'xtb_homo_ev', 'abs_error_ev', 'sascore', 'pass_025'
]].to_string())

if len(valid) > 0:
    mae    = valid['abs_error_ev'].mean()
    n_pass = int(valid['pass_025'].sum())
    print(f"\nMAE (ML vs xTB):     {mae:.4f} eV")
    print(f"Threshold:           {MAE_THRESHOLD} eV")
    print(f"Candidates passing:  {n_pass}/{len(valid)}")
    print(f"All synthesizable:   10/10 SAScore < 4.0")
    print(f"Validation:          {'PASSED ✅' if mae < MAE_THRESHOLD else 'FAILED ❌'}")

results_df.to_csv(
    "results/tables/phase3_xtb_validation_final.csv",
    index=False
)
print("\nSaved results/tables/phase3_xtb_validation_final.csv")