import pandas as pd
import numpy as np
from econml.dml import LinearDML
from sklearn.ensemble import RandomForestRegressor

print("=== EconML (Python) vs DoubleML (R) Comparison ===\n")

# Load data
df = pd.read_csv("data/processed/master_acceptor_dataset.csv")
df = df[df['ewg_count'] > 0].copy()
df['meas_num'] = (df['measurement_type'] == 'experiment').astype(int)
print(f"Rows: {len(df)}")

# Define variables — same as R
confounders = ['mol_weight', 'n_arom_rings', 'conj_length',
               'halogen_count', 'meas_num']
treatment = 'ewg_weighted'

outcomes = {
    'homo_ev':    'HOMO',
    'lumo_ev':    'LUMO',
    'bandgap_ev': 'Bandgap'
}

# R DoubleML results for comparison
r_results = {
    'homo_ev':    -0.01201265,
    'lumo_ev':    -0.02480654,
    'bandgap_ev': -0.01784393
}

results = []

for col, name in outcomes.items():
    Y = df[col].values
    T = df[treatment].values
    W = df[confounders].values

    model = LinearDML(
        model_y = RandomForestRegressor(n_estimators=500, min_samples_leaf=5,
                                        random_state=42),
        model_t = RandomForestRegressor(n_estimators=500, min_samples_leaf=5,
                                        random_state=42),
        cv      = 5,
        random_state = 42
    )
    model.fit(Y, T, X=None, W=W)

    ate   = model.ate()
    ci    = model.ate_interval(alpha=0.05)
    r_ate = r_results[col]
    diff  = abs(ate - r_ate)
    agree = diff < 0.005

    print(f"\n=== EWG -> {name} ===")
    print(f"EconML ATE:  {ate:.6f} eV")
    print(f"DoubleML ATE:{r_ate:.6f} eV")
    print(f"Difference:  {diff:.6f} eV")
    print(f"Within 0.005 eV threshold: {agree}")

    results.append({
        'outcome':     col,
        'econml_ate':  ate,
        'doubleml_ate': r_ate,
        'difference':  diff,
        'within_0005': agree,
        'ci_lower':    ci[0],
        'ci_upper':    ci[1]
    })

results_df = pd.DataFrame(results)
results_df.to_csv(
    "results/tables/econml_vs_doubleml_comparison.csv",
    index=False
)
print("\n=== Summary ===")
print(results_df[['outcome','econml_ate','doubleml_ate','difference','within_0005']])
print("\nSaved results/tables/econml_vs_doubleml_comparison.csv")