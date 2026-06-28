# Causal ML for Organic Photovoltaic Molecular Design

**Author:** Samuel Bizimana | JUNIA ISEN  
**Supervisor:** Prof Kekel N'KNOUN  
**Target Journal:** Nature Machine Intelligence

## What this project does
Applies causal inference methods to 12,138 OPV molecules to identify 
what *causes* solar cell efficiency — not just what correlates with it.

## Key Finding
A previously unknown electronic boundary at HOMO ≈ −5.5 eV where 
bandgap engineering transitions from beneficial to detrimental for 
power conversion efficiency (PCE). Confirmed statistically (p = 0.010).

## Methods Used
- Double Machine Learning (Chernozhukov et al., 2018)
- Sensitivity Analysis — robustness value = 0.469 (Cinelli & Hazlett, 2020)
- Causal Forests, 4000 trees (Wager & Athey, 2018)

## Progress
- Phase 0  Reproducible environment (R 4.6.0 + Python 3.10)
- Phase 1  Dataset — 12,138 OPV molecules
- Phase 2  Causal inference + CATE maps

Next: 
- Phase 3 🔄 Counterfactual molecular design (in progress)
- Phase 4 ⏳ Manuscript submission

## Reproduce
```r
renv::restore()  # R packages
```
```bash
conda env create -f environment.yml  # Python packages
```
