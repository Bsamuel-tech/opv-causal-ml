
## Data Availability

### Experimental dataset (599 molecules)
- Source: Published OPV literature, manually curated donor-acceptor pairs
- Original file: Data_Merged_with_SMILES.xlsx (1,573 donor-acceptor pairs)
- Processing: Acceptor molecules extracted, deduplicated by canonical SMILES
- Properties: HOMO, LUMO, optical bandgap measured experimentally
- Filter applied: ewg_count > 0 for causal analysis (598 molecules used)

### Harvard Clean Energy Project (14,930 molecules)
- Source: https://www.molecularspace.org
- Reference: Hachmann et al. (2011), J. Phys. Chem. Lett., 2(17), 2241-2251
- Access: Full database downloaded as moldata.csv (563 MB, not tracked in Git)
- Filters applied:
  - Random sample of 15,000 from 2,322,849 total molecules (seed=42)
  - HOMO: -8 to -3 eV
  - LUMO: -6 to -1 eV
  - Bandgap: 0.5 to 4.0 eV
  - Acceptor structural filter: LUMO < -2.0 eV, Bandgap < 3.5 eV, HOMO < -4.5 eV
  - 14,828/14,930 sampled molecules pass acceptor filter (99.32%)
- Properties: DFT-computed at B3LYP level

### Combined dataset
- Final: 15,529 molecules (599 experimental + 14,930 CEP)
- Saved: data/processed/master_acceptor_dataset.csv
- measurement_type column distinguishes experimental from DFT_B3LYP


## Model Files

results/models/causal_forest_ewg_homo.rds — trained causal forest (4,000 trees, 5MB)
This file is committed to GitHub and can be loaded directly.

causal_forest.rds (115MB) was a duplicate model file generated during development.
It has been deleted. To retrain from scratch run scripts/phase2_causal/causal_forest.R.
Training takes approximately 2 minutes on a standard laptop.
