# Contributing

Thanks for your interest in PRISM.

## How to report issues
- Please open a GitHub Issue with:
  - your OS / Python version
  - the command you ran
  - the full error message (copy/paste)
  - minimal input files if possible (avoid large trajectories)

## Pull requests
- Keep changes focused and well-described.
- Please do not commit large artifacts (MD trajectories, container images, large zip files).
- If you add dependencies, update:
  - conda/environment.yml
  - requirements.txt
  - apptainer/PRISM.def (if relevant)

## Scope note
- Step 6 production MD runs typically require HPC resources and are not expected to run in GitHub CI.
- The repository aims to keep analysis/post-processing scripts reproducible and testable.
