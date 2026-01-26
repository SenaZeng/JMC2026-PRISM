core/ – Workflow Scripts (Step1–Step6)

Step 1 – Fragment-based generation (FRATTVAE)
- Script: step1_vae.py
- Input: pretrained model files under results/pretrained_model/ (not included)
- Output: generated SMILES list (CSV)

Step 2 – Surrogate activity prediction (Uni-Mol + RF regressor)
- Script: step2_surrogate.py
- Input: curated Mpro bioactivity CSV (public; see data/README.md)
- Output: activity_predictor.pkl and prediction CSV

Step 2b – hERG liability classifier (auxiliary ADMET filter)
- Script: step2b_train_herg_model.py
- Input: herg_tdc_full.csv (public; see data/README.md)
- Output: herg_rf_model.pkl

Step 4 – ADMET + physics re-ranking
- Step4a: ADMET gating (+ hERG inference)
- Step4b: PySCF refinement
- Step4c: merge utilities

Step 5 – Broad-spectrum docking (worst-case score)
- Step5a: docking
- Step5b: merge utilities

Step 6 – MD preparation, execution, and analysis
- Step6c/6d: MD automation (HPC scheduler-specific options)
- Step6e: trajectory analysis + composite score
