
## Hardcoded Path Notice

Scripts in scripts/phase1_data/ contain hardcoded Windows paths:
  C:/Users/Samuel Bizimana/OneDrive/Desktop/Research Training/

To run on a different machine:
1. Clone the repository
2. Set your working directory to the repository root
3. Run scripts from the repository root (not from within subfolders)
4. All data/ and results/ paths are relative to the repository root

The corrected pipeline scripts (step1, step2, fix1, fix2) use relative
paths and will work from any machine when run from the repository root.
The hardcoded paths in step1_build_dataset.py and step2_expand_dataset.py
will be refactored in the next iteration to use os.path or a config file.

