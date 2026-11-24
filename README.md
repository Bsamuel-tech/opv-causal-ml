# Causal Machine Learning for Organic Acceptor Design

**Author:** Samuel Bizimana  
**Institution:** JUNIA ISEN  
**Supervisor:** Dr. Kekeli N'KONOU  
**Target Journal:** Nature Machine Intelligence  
**GitHub Repository:** https://github.com/Bsamuel-tech/opv-causal-ml  
**Pipeline Documentation:** https://github.com/Bsamuel-tech/opv-causal-ml/blob/main/PIPELINE.md

## 1. Project Overview

This project investigates the relationship between electron-withdrawing groups (EWGs) and the electronic properties of organic acceptor molecules using causal machine learning.

The assembled dataset contains **15,403 molecular records**, while the primary causal analysis is restricted to **598 experimental molecules** containing at least one EWG. The treatment variable is a continuous Hammett-weighted EWG score, used to investigate how EWG substitution relates to frontier orbital energies and molecular bandgaps.

The workflow covers dataset construction, molecular feature engineering, causal estimation, sensitivity analysis, heterogeneous treatment-effect analysis, counterfactual molecular design, and quantum chemistry validation.

The project is organised as a reproducible pipeline, with analysis scripts, processed data, result tables, figures, and model outputs maintained in the repository.

## 2. Key Results

### 2.1 Double Machine Learning estimates

Double Machine Learning (DML) estimates suggest a negative relationship between the Hammett-weighted EWG score and both HOMO and LUMO energies.

| Outcome | DML estimate (eV) | IID p-value | Cluster-robust p-value |
|---|---:|---:|---:|
| HOMO | −0.01828 | 0.042 | 0.052 |
| LUMO | −0.01633 | 0.049 | 0.051 |
| Bandgap | −0.00228 | 0.814 | 0.841 |

For HOMO, the IID 95% confidence interval is [−0.03589, −0.00066], while the scaffold cluster-robust interval is [−0.03671, +0.00016].

For LUMO, the corresponding intervals are [−0.03260, −0.00006] and [−0.03273, +0.00007].

The HOMO and LUMO estimates therefore have consistent negative directions, but their cluster-robust confidence intervals marginally cross zero. Neither outcome provides statistically significant evidence at the 5% level under cluster-robust inference. After Bonferroni correction for three outcomes, the adjusted p-values are approximately 0.13 for HOMO and 0.15 for LUMO. No statistically significant effect is detected for the bandgap.

The estimates were obtained using scaffold-grouped five-fold cross-fitting with cyclic round-robin fold assignment, across 20 repetitions. The analysis included 598 experimental molecules distributed across 288 unique molecular scaffolds. The reported results are from the committed analysis outputs and have not been independently reproduced.

Post-hoc scaffold cluster-robust standard errors were calculated using a sandwich estimator, following Cameron and Miller (2015), with row-wise averaging across the 20 repetitions. The results are available in `results/tables/dml_cluster_robust.csv`.

### 2.2 Instrumental variable analysis

An instrumental variable (IV) analysis was conducted using `nonarom_cc` as the sole instrument.

The first-stage F-statistic was 215.448. The Wu–Hausman tests returned p-values of 0.036 for HOMO and 0.026 for LUMO, indicating statistically detectable differences between the OLS and IV estimates under the assumptions of the test.

However, these results do not establish that the instrument is valid. The exclusion restriction remains chemically uncertain because non-aromatic carbon–carbon double bonds may directly influence molecular conjugation and frontier orbital energies.

The IV coefficients, approximately −0.066 eV, are around 3.7 times larger than the corresponding DML estimates. Although their directions are consistent, the difference in magnitude and uncertainty surrounding instrument validity prevent treating the IV results as independent confirmation of the DML findings.

The IV analysis is therefore reported as supplementary evidence in Section S3. DML remains the primary estimation approach.

### 2.3 Confounding and scaffold leakage

`halogen_count` is retained as a confounder because it captures structural information beyond its weighted contribution to the Hammett EWG score.

A bootstrap leakage analysis was conducted to assess the effect of molecular scaffold overlap on predictive performance.

Across the full dataset, the mean leakage gap was 0.0194 over 20 random seeds, with all observed gaps below the predefined 0.15 threshold.

For the 598-molecule causal subset, the mean gap was substantially larger at 0.316, with a 95% confidence interval of [0.164, 0.486].

This difference indicates that random splitting can overestimate predictive performance when related molecular scaffolds appear in both training and test sets. It supports the use of scaffold-grouped cross-fitting for the primary causal analysis.

### 2.4 Nuisance-model diagnostics

The nuisance-model diagnostics were:

- Treatment model (`ml_l`): R² = 0.181; RMSE = 0.133 eV.
- Outcome model (`ml_m`): R² = 0.735.

These metrics describe the predictive performance of the nuisance models used in the DML procedure.

### 2.5 Causal forest analysis

Causal forests were used to estimate conditional average treatment effects (CATEs) and explore potential variation in the effect of the EWG score across molecular structures.

The calibration test for the mean forest prediction was statistically significant (one-sided p = 0.0039). However, the differential prediction test did not provide statistically significant evidence of treatment-effect heterogeneity (one-sided p = 0.629).

The results therefore do not establish meaningful variation in treatment effects across the analysed molecules. The CATE maps are presented as exploratory descriptions of the estimated effects rather than evidence of confirmed molecular subgroups.

### 2.6 Sensitivity analysis

Sensitivity analysis was performed using `sensemakr` in R, with molecular weight as the benchmark confounder.

The reported HOMO sensitivity statistics were:

- t-value: −5.34
- Partial R²: 0.046
- Robustness value (point estimate): 19.7%
- Robustness value (at α = 0.05): 12.95%

These results describe the sensitivity of the corresponding OLS estimate to unobserved confounding. They should not be interpreted as a sensitivity analysis of the DML estimator.

### 2.7 Counterfactual design and quantum chemistry validation

The Phase 3 design workflow uses molecular features and estimated treatment effects to explore candidate modifications to organic acceptor molecules.

The resulting candidates were evaluated using xTB 6.6.1 with GFN2-xTB calculations and a fixed conformer seed of 42.

The validation results were:

- ML-predicted HOMO shift MAE: 0.15413 eV.
- Zero-shift baseline MAE: 0.15362 eV.
- Candidates within the predefined 0.25 eV shift-error threshold: 8/10 for both approaches.
- Candidates passing the SAScore filter (SAScore < 4.0): 10/10.

The ML model did not outperform the zero-shift baseline on this ten-candidate validation set. Consequently, the current results do not demonstrate improved predictive performance over simply assuming no change in HOMO energy.

The design workflow should therefore be considered exploratory. The validation results provide a starting point for further development rather than evidence of a validated molecular design system.

## 3. Analysis Scope and Restrictions

The primary causal analysis applies two filters:

```python
measurement_type == "experiment"
ewg_count > 0
```

This restricts causal estimation to **598 experimental molecules** containing at least one EWG group.

Nearly all CEP molecules have zero EWGs under the project's SMARTS-based definition. These molecules contribute to the broader dataset and scaffold diversity but are not included in the primary causal estimation.

The estimand is therefore the average treatment effect (ATE) among EWG-containing experimental organic acceptors. It does not represent the full population of 15,403 molecular records.

Further details on dataset composition and filtering are provided in `results/tables/supplementary_information.csv`, Section S0.

## 4. Methods and Implementation

The project combines cheminformatics, causal inference, machine learning, and quantum chemistry.

| Component | Method or tool |
|---|---|
| Molecular processing | RDKit (Python) |
| EWG identification | SMARTS-based molecular pattern matching |
| Causal estimation | Double Machine Learning using DoubleML in R |
| Cross-fitting | Scaffold-grouped five-fold cyclic round-robin assignment |
| Cluster-robust inference | Sandwich estimator following Cameron and Miller (2015) |
| Instrumental variable analysis | IV estimation using `nonarom_cc` |
| Sensitivity analysis | `sensemakr` in R |
| Causal forest analysis | `grf` in R, 4,000 trees and 288 scaffold clusters |
| Python causal comparison | EconML |
| Counterfactual design | Molecular modification and candidate generation |
| Quantum chemistry validation | xTB 6.6.1, GFN2-xTB |

DML uses random forest nuisance learners with five-fold scaffold-grouped cross-fitting and 20 repetitions. Cyclic round-robin fold assignment ensures that all folds receive observations.

Scaffold cluster-robust standard errors are calculated after estimation using a sandwich estimator. The IV analysis is supplementary because the chemical validity of the exclusion restriction remains uncertain.

Quantum chemistry calculations use a fixed conformer seed to improve reproducibility.

## 5. Repository Structure

The repository is organised into four main stages.

### Phase 1: Dataset construction

`scripts/phase1_data/`

Contains:

- `step1_build_dataset.py`
- `step2_expand_dataset.py`
- `fix1_add_columns.py`
- `fix2_scaffold_split.py`
- `fix3_acceptor_filter.py`
- `fix4_bootstrap_leakage.py`
- `ewg_utils.py` — shared EWG processing module.

### Phase 2: Causal analysis

`scripts/phase2_causal/`

Contains:

- `dml_analysis.R`
- `iv_analysis.R`
- `sensitivity_analysis.R`
- `causal_forest.R`
- `econml_comparison.py`

### Phase 3: Counterfactual design and validation

`scripts/phase3_design/`

Contains:

- `step1_counterfactual_design.py`
- `step2_xtb_validation.py`

The xTB validation script uses a fixed conformer seed for reproducibility.

### Results and outputs

`results/figures/`

Contains publication figures at 300 dpi, including:

- `cate_ewg_homo_distribution.png`
- `cate_ewg_vs_molweight.png`

`results/tables/`

Contains the analysis outputs and supplementary information, including:

- `dml_ewg_corrected.csv`
- `dml_fold_diagnostics.csv`
- `dml_cluster_robust.csv`

`results/models/`

Contains the trained causal forest model:

- `causal_forest_ewg_homo.rds` (approximately 5 MB)

### Processed data

`data/processed/`

Contains:

- `master_acceptor_dataset.csv` — the assembled dataset of 15,403 molecular records.

GitHub's browser preview may not display files larger than 1 MB. To access these files, open the file on GitHub and select **View raw**, or clone the repository locally.

## 6. Project Status

<details>
<summary>Phase 0 — Environment and reproducibility</summary>

**Status: Complete**

- R 4.6.0 and Python 3.10.
- R dependencies recorded in `renv.lock`.
- Python dependencies recorded in `requirements.txt` and `environment.yml`.

</details>

<details>
<summary>Phase 1 — Dataset construction</summary>

**Status: Complete**

- 15,403 assembled molecular records.
- 599 experimental acceptors after dataset construction and filtering.
- CEP acceptor filter pass rate: 99.32%.
- Full-dataset bootstrap leakage gap mean: 0.0194 across 20 seeds.
- Causal subset leakage analysis completed separately.

</details>

<details>
<summary>Phase 2 — Causal analysis</summary>

**Status: Complete**

- DML estimates obtained from committed code; not independently reproduced.
- Analysis run: 19 October 2026.
- 598 experimental molecules and 288 scaffolds.
- Five-fold scaffold-grouped cross-fitting, 20 repetitions.
- HOMO: IID p = 0.042; cluster-robust p = 0.052.
- LUMO: IID p = 0.049; cluster-robust p = 0.051.
- Bandgap: no statistically significant effect.
- Bonferroni-adjusted HOMO and LUMO p-values: approximately 0.13 and 0.15.
- IV analysis and instrument-validity limitation documented.
- Sensitivity analysis and CATE maps completed.

</details>

<details>
<summary>Phase 3 — Counterfactual design and validation</summary>

**Status: Complete — exploratory**

- Counterfactual molecular design workflow implemented.
- xTB validation completed for ten candidates.
- ML shift-prediction MAE: 0.15413 eV.
- Zero-shift baseline MAE: 0.15362 eV.
- Eight of ten candidates pass the 0.25 eV shift-error threshold under both approaches.
- All ten retained candidates pass the SAScore < 4.0 filter.
- ML did not outperform the zero-shift baseline on this validation set.

</details>

<details>
<summary>Phase 4 — Figures and supplementary material</summary>

**Status: Complete**

- Four publication figures.
- Seven-section supplementary information.
- Analysis tables and supporting outputs included in the repository.

</details>

## 7. Data Availability and Reproducibility

Dataset provenance, source information, processing decisions, and data availability details are documented in:

- `DATA_AVAILABILITY.md`
- `PATHS_NOTICE.md`

The pipeline documentation provides additional details about the processing and analysis stages:

https://github.com/Bsamuel-tech/opv-causal-ml/blob/main/PIPELINE.md

## 8. Reproducing the Analysis

### Restore the R environment

```r
renv::restore()
```

### Run the primary DML analysis

```r
source("scripts/phase2_causal/dml_analysis.R")
```

### Create the Python environment

```bash
conda env create -f environment.yml
conda activate causal-mol
```

### Run the bootstrap leakage analysis

```bash
python scripts/phase1_data/fix4_bootstrap_leakage.py
```

For the complete workflow, follow the instructions in `PIPELINE.md`.

---

## Acknowledgements

This work was conducted at JUNIA ISEN under the supervision of Dr. Kekeli N'KONOU.
