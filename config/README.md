# Configuration

PRISM uses a simple YAML configuration to store site-specific paths and external tool locations.

## Quick start
1) Copy the template:
   cp config/config.template.yaml config/config.yaml
2) Edit config/config.yaml to set:
   - paths.workdir (recommended outside the repo for large outputs)
   - optional: paths.models_dir / paths.datasets_dir
   - optional: executables.gromacs_gmx / xtb / docking engine paths

## Notes
- Do NOT commit config/config.yaml if it contains local paths, credentials, or site-specific settings.
- The repository ships config/config.template.yaml as a safe starting point.
