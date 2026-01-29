#!/usr/bin/env bash
set -euo pipefail

# Build PRISM container image from definition file
# Usage:
#   bash apptainer/build.sh
# Output:
#   prism.sif (in repo root)

apptainer build prism.sif apptainer/PRISM.def
