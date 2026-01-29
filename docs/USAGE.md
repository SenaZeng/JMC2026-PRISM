# Usage guide (minimal)

This repository provides modular PRISM scripts under core/. Many modules require external tools, models, or site-specific paths (see config/).

## Quick verification (recommended)
- Lightweight: create conda env and run:
  python examples/first_test.py
  (syntax-only; does not execute heavy workflows)

## Running modules
1) Prepare local configuration
- Copy config template:
  cp config/config.template.yaml config/config.yaml
- Edit config/config.yaml to define workdir and external executables/paths if needed.

2) Typical execution scope
- Steps 1–5: intended to run in a containerized environment (Apptainer/Singularity recipe provided).
- Step 6 production MD: typically requires HPC/site-installed GROMACS.
- Step 6 analysis: post-processing and scoring scripts can be applied to trajectories produced externally.

## Script map (core/)
- Step 1: step1_vae.py
- Step 2: step2_surrogate.py, step2b_train_herg_model.py
- Step 3: step3a_optimizer.py, step3b_run_dft.py, step3c_dft_refine.py
- Step 4: step4a_admet.py, step4b_final_pyscf.py
- Step 5: step5a_docking.py, step5b_utils_merge.py
- Step 6: step6a_pretop10.py, step6c_run_gromacs.py, step6d_run_gromacs_run.py, step6e_analysis.py

## Notes
- Some scripts may require large inputs (models, receptor files, trajectories) that are intentionally not committed (see MANIFEST.md).
- For container usage:
  bash apptainer/build.sh
  bash apptainer/run.sh python3 core/<script>.py --help
