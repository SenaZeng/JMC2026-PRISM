# Minimal runnable examples (PRISM)

This folder provides a lightweight first test to verify that the repository is callable after cloning.
It does not execute the full PRISM workflow or require large datasets/HPC.

## first test (Conda)
1) Create env:
   conda env create -f conda/environment.yml
2) Activate:
   conda activate prism
3) Run:
   python examples/first_test.py

Expected outcome: the script checks core/step1.py ... core/step6.py and attempts to run each with '--help' (or falls back to no-arg execution).
