PRISM (JMC2026) – Code and Environment for Broad-Spectrum Coronavirus Mpro Inhibitor Design

This repository provides representative code and environment specifications supporting the PRISM workflow described in our Journal of Medicinal Chemistry (JMC) submission.

Scope
- Fragment-based molecular generation (FRATTVAE wrapper)
- Surrogate-guided activity prediction (Uni-Mol representation + RandomForest regressor)
- Auxiliary safety filtering (hERG liability classifier; RandomForest + Morgan fingerprints)
- Physics-informed re-ranking (xTB/DFT descriptors; PySCF refinement)
- Multi-target docking (SARS-CoV-2 / SARS-CoV-1 / MERS Mpro) with worst-case scoring
- Molecular dynamics automation and analysis (QC gate, RMSD/Rg/mindist/H-bonds, composite robustness score)

Repository Structure
- core/               Core workflow scripts (Step1–Step6)
- env/                Environment files (conda/pip)
- data/README.md       Data provenance and how to obtain public datasets
- results/README.md    Output file schema and example column definitions
- docs/               Optional: additional notes for reproducibility

Quick Start (minimal)
1) Create environment
   - conda env create -f env/environment.yml
   - conda activate [ENV_NAME]

2) Run representative workflow steps (examples; adjust paths as needed)
   - Step1: generation
     python core/step1_vae.py
   - Step2: surrogate model (train or load)
     python core/step2_surrogate.py
   - Step2b: hERG model (optional; trains from public dataset)
     python core/step2b_train_herg_model.py
   - Step4–6: see core/README.md for step-specific commands and I/O

Reproducibility Notes
- HPC-specific job scheduling scripts and large intermediate artifacts (e.g., full MD trajectories) are not included.
- All algorithmic steps, parameter choices, and decision criteria are documented in the manuscript Methods and Supporting Information.

Data Availability
- Bioactivity datasets for surrogate and hERG models are from public sources (see data/README.md).
- Protein structures were retrieved from the Protein Data Bank (PDB). PDB IDs used: [LIST_PDBS].

License
[MIT / Apache-2.0 / TBD]

Contact
For questions related to reproduction of the workflow, please contact: [YOUR_EMAIL or omit]

