import pandas as pd
import numpy as np

val = pd.read_csv("results/tables/phase3_xtb_validation_final.csv")
print("Validation molecules:", len(val))
print("Columns:", list(val.columns))

# The 8 calibration molecules are those used to compute the 4.6197 eV offset
# We bootstrap these to quantify calibration uncertainty
cal = val[['original_homo_ev', 'xtb_homo_ev']].dropna().copy()
print("Calibration pairs available:", len(cal))

# Raw offsets: experimental - raw_xtb
# raw_xtb = xtb_homo_ev - 4.6197 (reverse the calibration)
cal['raw_xtb'] = cal['xtb_homo_ev'] - 4.6197
cal['offset']  = cal['original_homo_ev'] - cal['raw_xtb']

print("\n=== Calibration offset per molecule ===")
print(cal[['original_homo_ev', 'raw_xtb', 'offset']].round(4))

print(f"\nMean offset: {cal['offset'].mean():.4f} eV")
print(f"Std offset:  {cal['offset'].std():.4f} eV")

# Bootstrap the offset with 10000 resamples
np.random.seed(42)
n_boot = 10000
boot_offsets = []
for _ in range(n_boot):
    sample = cal['offset'].sample(n=len(cal), replace=True)
    boot_offsets.append(sample.mean())

boot_offsets = np.array(boot_offsets)
ci_lower = np.percentile(boot_offsets, 2.5)
ci_upper = np.percentile(boot_offsets, 97.5)

print(f"\n=== Bootstrap calibration uncertainty (n={n_boot} resamples) ===")
print(f"Bootstrap mean:  {boot_offsets.mean():.4f} eV")
print(f"Bootstrap std:   {boot_offsets.std():.4f} eV")
print(f"95% CI:          [{ci_lower:.4f}, {ci_upper:.4f}] eV")
print(f"CI width:        {ci_upper - ci_lower:.4f} eV")

# Propagate into MAE
mae_central = val['abs_error_ev'].mean()
print(f"\n=== MAE with bootstrap calibration uncertainty ===")
print(f"MAE central:     {mae_central:.4f} eV")
print(f"MAE lower bound: {(val['abs_error_ev'] - (ci_upper - boot_offsets.mean())).abs().mean():.4f} eV")
print(f"MAE upper bound: {(val['abs_error_ev'] + (ci_upper - boot_offsets.mean())).abs().mean():.4f} eV")
print(f"Threshold:       0.25 eV")

# Save
results = pd.DataFrame({
    'metric': [
        'N calibration molecules',
        'Calibration offset mean (eV)',
        'Calibration offset std (eV)',
        'Bootstrap 95% CI lower (eV)',
        'Bootstrap 95% CI upper (eV)',
        'Bootstrap CI width (eV)',
        'MAE central (eV)',
        'Threshold (eV)'
    ],
    'value': [
        len(cal),
        round(cal['offset'].mean(), 4),
        round(cal['offset'].std(), 4),
        round(ci_lower, 4),
        round(ci_upper, 4),
        round(ci_upper - ci_lower, 4),
        round(mae_central, 4),
        0.25
    ]
})
results.to_csv("results/tables/xtb_bootstrap_calibration.csv", index=False)
print("\nSaved results/tables/xtb_bootstrap_calibration.csv")
print(results.to_string(index=False))