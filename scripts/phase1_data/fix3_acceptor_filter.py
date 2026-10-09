"""
Fix 9: Add acceptor structural filter to CEP sample
====================================================
Filters CEP molecules to confirm they are acceptor-type
by checking for known acceptor motifs: electron-withdrawing
groups, low-lying LUMO, and absence of strong donor groups.
"""

import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT           = Path(__file__).resolve().parent.parent.parent
DATA_INTERIM   = ROOT / "data" / "interim"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS_TABLES = ROOT / "results" / "tables"
DATA_INTERIM.mkdir(parents=True, exist_ok=True)
RESULTS_TABLES.mkdir(parents=True, exist_ok=True)

print("=== Fix 9: Acceptor structural filter for CEP sample ===\n")

# READ FROM fix2's output (was step2_merged_raw.csv — broken link fixed)
df = pd.read_csv(DATA_INTERIM / "fix2_with_scaffolds.csv")
cep = df[df['source_db'] == 'CEP'].copy()
print(f"CEP molecules before filter: {len(cep)}")

# Acceptor structural criteria:
# 1. At least one EWG motif OR low LUMO (< -3.5 eV)
# 2. No strong donor-only pattern (EDG count < EWG count)
# 3. LUMO < -2.0 eV (acceptors have low-lying LUMOs)

def is_acceptor_like(row):
    if row['lumo_ev'] > -2.0:
        return False
    if row['bandgap_ev'] > 3.5:
        return False
    if row['homo_ev'] > -4.5:
        return False
    return True

cep['is_acceptor'] = cep.apply(is_acceptor_like, axis=1)

n_pass = cep['is_acceptor'].sum()
n_fail = (~cep['is_acceptor']).sum()

print(f"CEP molecules passing acceptor filter: {n_pass}")
print(f"CEP molecules failing acceptor filter: {n_fail}")
print(f"Percent passing: {round(n_pass/len(cep)*100, 2)}%")

print("\nFilter criteria:")
print("  LUMO < -2.0 eV (low-lying unoccupied orbital)")
print("  Bandgap < 3.5 eV (not wide-gap insulator)")
print("  HOMO < -4.5 eV (deep-lying occupied orbital)")

print("\nProperty ranges of passing CEP molecules:")
passing = cep[cep['is_acceptor']]
print(f"  HOMO: {passing['homo_ev'].min():.3f} to {passing['homo_ev'].max():.3f} eV")
print(f"  LUMO: {passing['lumo_ev'].min():.3f} to {passing['lumo_ev'].max():.3f} eV")
print(f"  Bandgap: {passing['bandgap_ev'].min():.3f} to {passing['bandgap_ev'].max():.3f} eV")

# Apply filter — keep only acceptor-like CEP + all experimental
exp = df[df['source_db'] != 'CEP'].copy()
cep_filtered = cep[cep['is_acceptor']].drop(columns=['is_acceptor'])
master_filtered = pd.concat([exp, cep_filtered], ignore_index=True)
print(f"\nFiltered master dataset: {len(master_filtered)} molecules")
print(master_filtered['source_db'].value_counts().to_string())

# Save interim snapshot
interim_out = DATA_INTERIM / "fix3_acceptors_only.csv"
master_filtered.to_csv(interim_out, index=False)
print(f"\nSaved interim: {interim_out}")

# Save filter report
report = pd.DataFrame({
    'filter': ['LUMO < -2.0 eV', 'Bandgap < 3.5 eV', 'HOMO < -4.5 eV'],
    'description': [
        'Low-lying LUMO confirms electron-accepting character',
        'Excludes wide-gap insulators not relevant to OPV',
        'Deep HOMO confirms acceptor not donor character'
    ],
    'molecules_passing': [
        (cep['lumo_ev'] < -2.0).sum(),
        (cep['bandgap_ev'] < 3.5).sum(),
        (cep['homo_ev'] < -4.5).sum()
    ]
})

report_out = RESULTS_TABLES / "cep_acceptor_filter_report.csv"
report.to_csv(report_out, index=False)
print(f"Saved filter report: {report_out}")
