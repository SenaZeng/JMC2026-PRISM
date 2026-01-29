# PRISM (JMC2026) — Containerized Workflow for Broad-Spectrum Coronavirus Mpro Inhibitor Discovery

PRISM is a modular AI + physics pipeline for de novo molecular design and multi-target screening of coronavirus main protease (Mpro) inhibitors. To support computational reproducibility, PRISM is distributed as a containerized workflow (Singularity/Apptainer) covering Steps 1–5, with consistent post-processing and scoring utilities for MD trajectories (Step 6 analysis).

## Repository structure
- core/                Step 1–Step 6 Python entry scripts (step1.py ... step6.py)
- apptainer/           (to be added) Apptainer/Singularity definition + build/run helpers
- conda/               (to be added) Lightweight environment for quick inspection/partial runs
- examples/            (to be added) Minimal runnable examples (small inputs + smoke test)
- docs/                (optional) Additional documentation

## Pipeline overview (Steps 1–6)
1) Structure–Semantic Generation (Steps 1–2): fragment-based molecular assembly and latent-space optimization.
2) Surrogate and Electronic Screening (Steps 3–4): Uni-Mol–guided prioritization and multi-fidelity electronic-structure screening.
3) Multi-Target Docking (Step 5): automated receptor preparation and broad-spectrum docking analysis (SARS-CoV-2 / SARS-CoV-1 / MERS Mpro).
4) Molecular Dynamics (Step 6): production MD runs executed in a site-installed GROMACS environment (HPC). Trajectory post-processing scripts, MDAnalysis workflows, and uncertainty-aware PRISM scoring are provided here to enable reproducible analysis when applied to resulting trajectories.

## Project status
- CI: GitHub Actions runs the minimal first test (syntax-only) on every push/PR.
- Citation: see CITATION.cff (Zenodo DOI will be added upon acceptance).
- License: MIT (see LICENSE).
- Releases: see CHANGELOG.md and RELEASE_CHECKLIST.md for versioned snapshots.


## Quick start (lightweight Conda; partial execution)
1) Create env:
   conda env create -f conda/environment.yml
2) Activate:
   conda activate prism
3) Run the first test:
   python examples/first_test.py

Note: this Conda environment is intended for rapid inspection and partial execution. Full dependency coverage for Steps 1–5 will be provided via the container recipe.

## Quick start (Apptainer/Singularity; full Steps 1–5 + MD analysis tools)
1) Build the container image:
   bash apptainer/build.sh
2) Run commands inside the container:
   bash apptainer/run.sh python3 core/step1.py --help

Example (run the first test inside the container):
   bash apptainer/run.sh python3 examples/first_test.py

## Reproducibility note
Steps 1–5 are designed to run inside the container for consistent dependencies. Step 6 production MD requires HPC resources and was executed using a site-installed GROMACS environment. However, the analysis layer (trajectory post-processing, MDAnalysis workflows, and PRISM scoring/robustness metrics) is included, so prioritization criteria are reproducible when the same analysis is applied to the resulting trajectories.

## Citation
Source code and container recipes: https://github.com/SenaZeng/JMC2026-PRISM
A versioned snapshot for the study will be archived on Zenodo (DOI to be assigned upon acceptance).
If you use PRISM, please cite this repository via the GitHub “Cite this repository” panel (CITATION.cff).

Versioned snapshot for this study:
- GitHub Release: v0.1.0

A Zenodo archive (DOI) will be added upon acceptance.

