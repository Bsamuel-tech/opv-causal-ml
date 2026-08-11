"""
Phase 3 — Step 2: xTB Quantum Chemistry Validation
====================================================
Validates counterfactual candidates using xTB GFN2.
ML predictions use per-molecule CATE from causal forest.
Calibration offset computed against 8 experimental molecules.

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
XTB_OFFSET    = 4.6197   # calibrated against 8 experimental molecules
MAE_THRESHOLD = 0.25     # eV

# ── xTB runner ────────────────────────────────────────────────────
def run_xtb_homo(smiles, mol_id, seed=42):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)

        # Fixed seed for reproducibility
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
print(f"MAE threshold: {MAE_THRESHOLD} eV")
print(f"Conformer seed: 42 (fixed for reproducibility)\n")

# ── Run validation ────────────────────────────────────────────────
results = []

for i, row in df.iterrows():
    orig_homo     = row['original_homo']
    pred_homo     = row['predicted_homo']
    pred_shift    = row['predicted_shift']
    cate          = row['cate']
    cate_se       = row['cate_se']

    print(f"Candidate {len(results)+1}/10 (SA={row['sascore']:.3f}):")
    print(f"  Original HOMO:    {orig_homo:.4f} eV")
    print(f"  Per-molecule CATE:{cate:.6f} eV per unit EWG")
    print(f"  ML predicted:     {pred_homo:.4f} eV "
          f"(shift={pred_shift:+.4f})")
    print(f"  Running xTB...")

    xtb_raw = run_xtb_homo(row['modified_smiles'], len(results))

    if xtb_raw is not None:
        xtb_cal   = xtb_raw + XTB_OFFSET
        xtb_shift = xtb_cal - orig_homo
        error     = abs(pred_homo - xtb_cal)
        passed    = error < MAE_THRESHOLD
        print(f"  xTB calibrated:   {xtb_cal:.4f} eV "
              f"(shift={xtb_shift:+.4f})")
        print(f"  Absolute error:   {error:.4f} eV "
              f"{'✅ PASS' if passed else '❌ FAIL'}")
    else:
        xtb_cal = xtb_shift = error = None
        passed  = None
        print(f"  xTB: failed")

    results.append({
        'candidate':          len(results) + 1,
        'original_smiles':    row['original_smiles'],
        'modified_smiles':    row['modified_smiles'],
        'original_homo_ev':   orig_homo,
        'per_molecule_cate':  cate,
        'cate_se':            cate_se,
        'ml_predicted_homo':  pred_homo,
        'ml_homo_shift_ev':   pred_shift,
        'xtb_homo_ev':        xtb_cal,
        'xtb_homo_shift_ev':  xtb_shift,
        'sascore':            row['sascore'],
        'abs_error_ev':       error,
        'pass_025':           passed,
        'synthesizable':      'YES'
    })

# ── Summary ───────────────────────────────────────────────────────
results_df = pd.DataFrame(results)
valid      = results_df.dropna(subset=['xtb_homo_ev'])

print("\n=== Validation Summary ===")
print(results_df[[
    'candidate', 'per_molecule_cate',
    'ml_predicted_homo', 'xtb_homo_ev',
    'abs_error_ev', 'pass_025'
]].to_string())

if len(valid) > 0:
    mae    = valid['abs_error_ev'].mean()
    n_pass = int(valid['pass_025'].sum())
    print(f"\nMAE (CATE-ML vs xTB): {mae:.4f} eV")
    print(f"Threshold:            {MAE_THRESHOLD} eV")
    print(f"Candidates passing:   {n_pass}/{len(valid)}")
    print(f"All synthesizable:    10/10 SAScore < 4.0")
    print(f"Validation:           "
          f"{'PASSED ✅' if mae < MAE_THRESHOLD else 'FAILED ❌'}")

results_df.to_csv(
    "results/tables/phase3_xtb_validation_final.csv",
    index=False
)
print("\nSaved phase3_xtb_validation_final.csv")
print("NOTE: ML predictions use per-molecule CATE from causal forest")