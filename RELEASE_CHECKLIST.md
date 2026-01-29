# Release checklist (GitHub Release + Zenodo)

## Before tagging
- [ ] CI is green on main.
- [ ] README Quick start commands are up to date.
- [ ] CITATION.cff author list and metadata are correct.
- [ ] Version number is decided (e.g., v1.0.0) and matches CHANGELOG.
- [ ] - [ ] CITATION.cff authors are real names (no placeholders) and ordered correctly.
- [ ] README Citation section points to the intended release tag (e.g., v0.1.0).


## Create GitHub Release
- [ ] Create a tag (e.g., v1.0.0) on main.
- [ ] Create a GitHub Release from that tag.
- [ ] Paste a short release note summarizing changes (from CHANGELOG).

## Zenodo archive + DOI
- [ ] Zenodo GitHub integration enabled for this repository.
- [ ] Confirm Zenodo created an archive for the release and minted a DOI.

## After DOI is minted
- [ ] Update CITATION.cff with the Zenodo DOI.
- [ ] Update README “Citation” section with the DOI.
- [ ] Update manuscript reproducibility statement with the DOI (replace placeholder).
