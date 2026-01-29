# Repository manifest (what is included / excluded)

## Included in this repository
- Source code for PRISM modules (core/).
- Lightweight environments (conda/environment.yml, requirements.txt).
- Minimal first test and examples (examples/).
- Apptainer container recipe and helper scripts (apptainer/PRISM.def, build.sh, run.sh).
- Documentation and project health files (README, LICENSE, CITATION.cff, SECURITY, CONTRIBUTING, CHANGELOG).

## Intentionally excluded (too large / site-specific / regenerate)
- Production MD trajectories and large simulation outputs (e.g., *.xtc, *.trr, *.dcd).
- Prebuilt container images (*.sif).
- Large intermediate/result archives (e.g., *.zip, *.tar.gz).
- Site-installed HPC software stacks (e.g., GROMACS production runs).

## Notes
- Steps 1–5 are designed to run inside the container for consistent dependencies.
- Step 6 production MD is typically executed in an HPC environment (site-installed GROMACS).
- Trajectory post-processing scripts and PRISM scoring logic are included to enable reproducible analysis when applied to resulting trajectories.
- A versioned snapshot of the repository will be archived on Zenodo upon acceptance (DOI to be assigned).
