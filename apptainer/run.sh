#!/usr/bin/env bash
set -euo pipefail

# Run commands inside the PRISM container image.
# Usage examples:
#   bash apptainer/run.sh python3 core/step1.py --help
#   bash apptainer/run.sh python3 examples/first_test.py

IMAGE="prism.sif"

if [[ ! -f "${IMAGE}" ]]; then
  echo "[ERROR] Container image not found: ${IMAGE}"
  echo "Build it first: bash apptainer/build.sh"
  exit 1
fi

apptainer exec "${IMAGE}" "$@"
