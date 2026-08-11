## Progress

-  Phase 0 : Reproducible environment (R 4.6.0 + Python 3.10)
-  Phase 1 : Dataset: 15,529 acceptor molecules, scaffold split confirmed
-  Phase 2 : DML + IV + sensitivity analysis + CATE maps
-  Phase 3 : Counterfactual design engine + xTB validation (MAE=0.220 eV)
-  Phase 4 : 4 publication figures + 6-section supplementary information

## Large files (kept locally only)

The following files exceed GitHub limits and are not tracked:

- `data/corrected_pipeline/master_acceptor_dataset.csv` — full dataset (563 MB)
- `data/raw/moldata.csv` — Harvard CEP raw data (563 MB)

To reproduce the dataset, run `scripts/phase1_data/step1_build_dataset.py`
followed by `scripts/phase1_data/step2_expand_dataset.py`.

## Reproduce

```r
renv::restore()  # restores all R packages
```

```bash
conda activate causal-mol  # Python environment
```

## References

- Chernozhukov et al. (2018). Double/debiased machine learning. Econometrics Journal.
- Wager & Athey (2018). Causal forests. JASA.
- Ertl & Schuffenhauer (2009). SAScore. Journal of Cheminformatics.
- Bannwarth et al. (2020). xTB methods. WIREs Computational Molecular Science.
- Cinelli & Hazlett (2020). Sensemakr. Journal of the Royal Statistical Society B.
- Hachmann et al. (2011). Harvard CEP. Journal of Physical Chemistry Letters.