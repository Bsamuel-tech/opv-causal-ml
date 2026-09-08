import pandas as pd
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.model_selection import cross_val_predict, KFold
from sklearn.preprocessing import StandardScaler

df = pd.read_csv("data/processed/master_acceptor_dataset.csv")
df_exp = df[df['source_db'] == 'your_data'].copy()
print("Experimental molecules available:", len(df_exp))

# Load xTB validation results
val = pd.read_csv("results/tables/phase3_xtb_validation_final.csv")
val = val.dropna(subset=['xtb_homo_ev', 'original_homo_ev'])
print("Molecules with both experimental and xTB HOMO:", len(val))

# Compute offset distribution across all available pairs
offsets = val['original_homo_ev'] - (val['xtb_homo_ev'] - 4.6197)
print("\n=== Calibration offset analysis ===")
print(f"Mean offset:   {offsets.mean():.4f} eV")
print(f"Std offset:    {offsets.std():.4f} eV")
print(f"95% CI:        [{offsets.mean() - 1.96*offsets.std():.4f}, "
      f"{offsets.mean() + 1.96*offsets.std():.4f}] eV")
print(f"Calibration offset used: 4.6197 eV")
print(f"Offset std as fraction of MAE threshold (0.25 eV): "
      f"{offsets.std()/0.25:.3f}")

# Propagate calibration uncertainty into MAE
val['xtb_homo_corrected'] = val['xtb_homo_ev']
val['abs_error'] = abs(val['ml_predicted_homo'] - val['xtb_homo_corrected'])
val['lower_error'] = abs(val['ml_predicted_homo'] -
                        (val['xtb_homo_corrected'] - 1.96*offsets.std()))
val['upper_error'] = abs(val['ml_predicted_homo'] -
                        (val['xtb_homo_corrected'] + 1.96*offsets.std()))

mae_central = val['abs_error'].mean()
mae_lower   = val['lower_error'].mean()
mae_upper   = val['upper_error'].mean()

print(f"\n=== MAE with propagated calibration uncertainty ===")
print(f"MAE (central):  {mae_central:.4f} eV")
print(f"MAE (lower CI): {mae_lower:.4f} eV")
print(f"MAE (upper CI): {mae_upper:.4f} eV")
print(f"Threshold:      0.25 eV")

if mae_upper < 0.25:
    print("Result: MAE passes threshold even at upper CI bound")
elif mae_central < 0.25:
    print("Result: MAE passes at central estimate but not at upper CI bound")
else:
    print("Result: MAE does not pass threshold")

# Literature-based uncertainty for experimental HOMO measurement methods
exp_uncertainty = {
    'UPS (ultraviolet photoelectron spectroscopy)': 0.05,
    'CV (cyclic voltammetry)': 0.07,
    'Optical spectroscopy': 0.03
}

print("\n=== Literature experimental uncertainty ===")
for method, unc in exp_uncertainty.items():
    print(f"  {method}: {unc} eV")

mean_exp_unc = np.mean(list(exp_uncertainty.values()))
justified_threshold = max(mean_exp_unc, 4.6197 * offsets.std())
print(f"\nJustified threshold (max of exp uncertainty and calibration uncertainty):")
print(f"  {justified_threshold:.4f} eV")
print(f"  Adopted threshold (0.25 eV) is {'conservative' if 0.25 > justified_threshold else 'liberal'} "
      f"relative to this estimate")

# Save
results = pd.DataFrame({
    'metric': [
        'Calibration offset (eV)',
        'Offset std (eV)',
        'CI lower (eV)',
        'CI upper (eV)',
        'MAE central (eV)',
        'MAE lower CI (eV)',
        'MAE upper CI (eV)',
        'Threshold (eV)',
        'Mean experimental uncertainty (eV)',
        'Justified threshold (eV)'
    ],
    'value': [
        round(4.6197, 4),
        round(offsets.std(), 4),
        round(offsets.mean() - 1.96*offsets.std(), 4),
        round(offsets.mean() + 1.96*offsets.std(), 4),
        round(mae_central, 4),
        round(mae_lower, 4),
        round(mae_upper, 4),
        0.25,
        round(mean_exp_unc, 4),
        round(justified_threshold, 4)
    ]
})
results.to_csv("results/tables/xtb_calibration_uncertainty.csv", index=False)
print("\nSaved results/tables/xtb_calibration_uncertainty.csv")
print(results.to_string(index=False))