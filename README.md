# Causal ML for Organic Acceptor Molecular Design

**Author:** Samuel Bizimana | JUNIA ISEN  
**Supervisor:** Dr. Kekeli N'KONOU  
**Target Journal:** Nature Machine Intelligence

## What this project does

Applies causal inference to 15,529 organic acceptor molecules to identify
what *causally* drives frontier orbital energies — not just what correlates.
Treatment variable: Hammett-weighted electron-withdrawing group (EWG) score.
Outcomes: HOMO energy, LUMO energy, optical bandgap.

## Key findings (Phases 0–2)

- EWG score causally lowers LUMO energy by −0.025 eV (DML, p=0.0008)
- IV estimation confirms direction: F-statistic=304, Sargan p=0.182
- Scaffold split generalisation gap = 0.050 (minimal data leakage)
- Causal forest calibration confirmed (p=0.0004)

## Methods

- Double Machine Learning — DoubleML R package
- Instrumental Variable estimation — ivreg R package  
- Sensitivity analysis — sensemakr R package
- Causal forests and CATE maps — grf R package
- Molecular processing — RDKit (Python)
- Quantum chemistry validation — xTB 6.7.1 (Phase 3)

## Project structure

scripts/phase1_data/    — dataset construction and scaffold split

scripts/phase2_causal/  — DML, IV, sensitivity, causal forest

scripts/phase3_design/  — counterfactual generator (coming)

data/processed/         — master_acceptor_dataset.csv (15,529 molecules)

results/figures/        — all publication figures

results/tables/         — all results tables



## Progress

- Phase 0 Reproducible environment (R 4.6.0 + Python 3.10)
- Phase 1  Dataset — 15,529 acceptor molecules, scaffold split
- Phase 2  DML + IV + sensitivity + CATE maps



- Phase 3 🔄 Counterfactual molecular design (next)
- Phase 4 ⏳ Manuscript submission

## Reproduce

```r
renv::restore()  # restores all R packages
```

```bash
conda activate causal-mol  # Python environment
```
