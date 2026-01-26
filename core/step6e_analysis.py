#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step 6e: MD trajectory analysis

- Prefer fastmdanalysis (FastMDAnalysis) if available
- Fallback to MDAnalysis if fastmdanalysis is unavailable or fails (when engine=auto)

Input CSV must contain:
  - an ID column (default tries: ID, Name, smiles, MolID)
  - topology file path column (default: top_path)
  - trajectory file path column (default: traj_path)

Typical usage:
  conda activate md_env
  python core/step6a_MD.py \
    --input_csv results/md_input.csv \
    --out_dir results/step6a_md \
    --engine auto \
    --top_col top_path \
    --traj_col traj_path \
    --sel_ligand "resname LIG" \
    --stride 10
"""

import argparse
import json
import os
import sys
from typing import Dict, Any, List, Optional

import pandas as pd


def _safe_mkdir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def _pick_id_col(df: pd.DataFrame, preferred: List[str]) -> str:
    for c in preferred:
        if c in df.columns:
            return c
    raise ValueError(
        f"Cannot find an ID column. Tried: {preferred}. "
        f"Available: {list(df.columns)[:60]}..."
    )


def _require_col(df: pd.DataFrame, col: str) -> None:
    if col not in df.columns:
        raise ValueError(f"Missing required column '{col}'. Available: {list(df.columns)}")


def _file_exists(path: str) -> bool:
    return isinstance(path, str) and len(path) > 0 and os.path.exists(path)


def try_import_fastmdanalysis():
    """
    FastMDAnalysis project installs as package name 'fastmdanalysis' (lowercase).
    Some users may try 'fastMDanalysis' (camelcase) by mistake; we provide fallback.

    Return imported module object or None.
    """
    try:
        import fastmdanalysis  # type: ignore
        return fastmdanalysis
    except Exception:
        pass

    # backward/alternative (rare)
    try:
        import fastMDanalysis  # type: ignore
        return fastMDanalysis
    except Exception:
        return None


def analyze_with_fastmd(
    fastmd_module,
    top_path: str,
    traj_path: str,
    out_dir: str,
    cfg: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Best-effort wrapper for fastmdanalysis.
    Since fastmdanalysis APIs can vary, we try common entrypoints.
    If unknown, return failure and allow fallback to MDAnalysis (engine=auto).
    """
    _safe_mkdir(out_dir)
    summary: Dict[str, Any] = {
        "engine": "fastmdanalysis",
        "top_path": top_path,
        "traj_path": traj_path,
        "status": "init",
    }

    try:
        # Try: fastmdanalysis.run(...)
        if hasattr(fastmd_module, "run"):
            fastmd_module.run(top=top_path, traj=traj_path, out=out_dir, **cfg)
            summary["status"] = "ok"
            summary["note"] = "fastmdanalysis.run executed"
            return summary

        # Try: fastmdanalysis.Analysis(...).run()
        if hasattr(fastmd_module, "Analysis"):
            ana = fastmd_module.Analysis(top_path, traj_path, out_dir=out_dir, **cfg)
            if hasattr(ana, "run"):
                ana.run()
                summary["status"] = "ok"
                summary["note"] = "fastmdanalysis.Analysis(...).run() executed"
                return summary

        summary["status"] = "failed"
        summary["error"] = (
            "fastmdanalysis imported but no recognized API entrypoint found "
            "(expected run or Analysis). Please adapt analyze_with_fastmd() to your version."
        )
        return summary

    except Exception as e:
        summary["status"] = "failed"
        summary["error"] = f"fastmdanalysis exception: {e}"
        return summary


def analyze_with_mdanalysis(
    top_path: str,
    traj_path: str,
    out_dir: str,
    cfg: Dict[str, Any]
) -> Dict[str, Any]:
    """
    MDAnalysis-based fallback: computes basic metrics:
      - Protein backbone RMSD (aligned)
      - Protein radius of gyration (Rg)
      - Ligand-protein minimum distance (binding persistence proxy)
      - Hydrogen bond events total (best-effort; selections may need tuning)
    """
    _safe_mkdir(out_dir)
    summary: Dict[str, Any] = {
        "engine": "MDAnalysis",
        "top_path": top_path,
        "traj_path": traj_path,
        "status": "init",
    }

    try:
        import numpy as np
        import MDAnalysis as mda
        from MDAnalysis.analysis import align
        from MDAnalysis.analysis.rms import RMSD
        from MDAnalysis.analysis.distances import distance_array
        from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
    except Exception as e:
        summary["status"] = "failed"
        summary["error"] = f"MDAnalysis import error: {e}"
        return summary

    sel_protein = cfg.get("sel_protein", "protein")
    sel_bb = cfg.get("sel_backbone", "protein and backbone")
    sel_lig = cfg.get("sel_ligand", "resname LIG")
    stride = int(cfg.get("stride", 1))
    ref_frame = int(cfg.get("ref_frame", 0))

    hb_donors_sel = cfg.get("hb_donors_sel", "protein")
    hb_acceptors_sel = cfg.get("hb_acceptors_sel", "protein")

    try:
        u = mda.Universe(top_path, traj_path)
    except Exception as e:
        summary["status"] = "failed"
        summary["error"] = f"Universe load error: {e}"
        return summary

    # Align to backbone
    try:
        ref = mda.Universe(top_path, traj_path)
        ref.trajectory[ref_frame]
        aligner = align.AlignTraj(u, ref, select=sel_bb, in_memory=False, filename=None, verbose=False)
        aligner.run(step=stride)
        summary["aligned"] = True
        summary["align_select"] = sel_bb
    except Exception as e:
        summary["aligned"] = False
        summary["align_error"] = str(e)

    # Protein backbone RMSD
    try:
        R = RMSD(u, u, select=sel_bb, ref_frame=ref_frame)
        R.run(step=stride)
        arr = R.results.rmsd  # frame, time(ps), rmsd(A)
        rmsd_bb = arr[:, 2].astype(float)

        summary["protein_bb_rmsd_A_mean"] = float(np.mean(rmsd_bb))
        summary["protein_bb_rmsd_A_median"] = float(np.median(rmsd_bb))
        summary["protein_bb_rmsd_A_p90"] = float(np.percentile(rmsd_bb, 90))

        pd.DataFrame(arr, columns=["frame", "time_ps", "rmsd_A"]).to_csv(
            os.path.join(out_dir, "protein_bb_rmsd.csv"), index=False
        )
    except Exception as e:
        summary["protein_bb_rmsd_error"] = str(e)

    # Protein Rg
    try:
        prot = u.select_atoms(sel_protein)
        rgs = []
        times = []
        for ts in u.trajectory[::stride]:
            rgs.append(float(prot.radius_of_gyration()))
            times.append(float(getattr(ts, "time", ts.frame)))
        rgs_arr = np.asarray(rgs, dtype=float)

        summary["protein_rg_A_mean"] = float(np.mean(rgs_arr))
        summary["protein_rg_A_p90"] = float(np.percentile(rgs_arr, 90))

        pd.DataFrame({"time": times, "Rg_A": rgs_arr}).to_csv(
            os.path.join(out_dir, "protein_rg.csv"), index=False
        )
    except Exception as e:
        summary["protein_rg_error"] = str(e)

    # Ligand-protein minimum distance
    try:
        prot = u.select_atoms(sel_protein)
        lig = u.select_atoms(sel_lig)
        summary["ligand_sel_used"] = sel_lig

        if len(lig) == 0:
            summary["ligand_detected"] = False
        else:
            summary["ligand_detected"] = True
            mindists = []
            times = []
            for ts in u.trajectory[::stride]:
                d = distance_array(lig.positions, prot.positions)
                mindists.append(float(np.min(d)))
                times.append(float(getattr(ts, "time", ts.frame)))
            md_arr = np.asarray(mindists, dtype=float)

            summary["lig_prot_mindist_A_mean"] = float(np.mean(md_arr))
            summary["lig_prot_mindist_A_p10"] = float(np.percentile(md_arr, 10))
            summary["lig_prot_mindist_A_p90"] = float(np.percentile(md_arr, 90))

            pd.DataFrame({"time": times, "min_dist_A": md_arr}).to_csv(
                os.path.join(out_dir, "lig_prot_min_dist.csv"), index=False
            )
    except Exception as e:
        summary["lig_prot_mindist_error"] = str(e)

    # Hydrogen bond events total (best-effort)
        # Hydrogen bond events total (best-effort; API differs across MDAnalysis versions)
    try:
        hb_dist = float(cfg.get("hb_distance", 3.5))
        hb_ang = float(cfg.get("hb_angle", 150.0))

        try:
            h = HydrogenBondAnalysis(
                u,
                donors_sel=hb_donors_sel,
                acceptors_sel=hb_acceptors_sel,
                hydrogens_sel=None,
                distance=hb_dist,
                angle=hb_ang,
            )
        except TypeError:
            # older/newer API variants
            h = HydrogenBondAnalysis(
                u,
                donors_sel=hb_donors_sel,
                acceptors_sel=hb_acceptors_sel,
                hydrogens_sel=None,
                d_a_cutoff=hb_dist,
                d_h_a_angle_cutoff=hb_ang,
            )

        h.run(step=stride)
        hb_count = len(h.results.hbonds) if hasattr(h, "results") else None
        summary["hbond_events_total"] = int(hb_count) if hb_count is not None else None

        if hasattr(h, "results") and hasattr(h.results, "hbonds"):
            pd.DataFrame(h.results.hbonds).to_csv(os.path.join(out_dir, "hbonds_raw.csv"), index=False)
    except Exception as e:
        summary["hbond_error"] = str(e)

    summary["status"] = "ok"
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_csv", required=True, help="CSV listing systems to analyze")
    ap.add_argument("--out_dir", required=True, help="Output directory for step6 analysis")
    ap.add_argument("--id_col", default=None, help="ID column; auto-pick if not set")
    ap.add_argument("--top_col", default="top_path", help="Topology path column name")
    ap.add_argument("--traj_col", default="traj_path", help="Trajectory path column name")
    ap.add_argument("--engine", default="auto", choices=["auto", "fastmd", "mdanalysis"],
                    help="Analysis engine: auto prefers fastmdanalysis if available")
    ap.add_argument("--stride", type=int, default=1, help="Frame stride")
    ap.add_argument("--ref_frame", type=int, default=0, help="Reference frame for alignment/RMSD")
    ap.add_argument("--sel_backbone", default="protein and backbone", help="Alignment/RMSD selection")
    ap.add_argument("--sel_protein", default="protein", help="Protein selection")
    ap.add_argument("--sel_ligand", default="resname LIG",
                    help="Ligand selection (adjust to your system; e.g., 'resname UNL' or 'resname MOL')")
    ap.add_argument("--hb_donors_sel", default="protein", help="Hbond donors selection (MDAnalysis)")
    ap.add_argument("--hb_acceptors_sel", default="protein", help="Hbond acceptors selection (MDAnalysis)")
    ap.add_argument("--hb_distance", type=float, default=3.5, help="Hbond distance cutoff (A)")
    ap.add_argument("--hb_angle", type=float, default=150.0, help="Hbond angle cutoff (deg)")
    ap.add_argument("--config_json", default=None, help="Optional JSON file for extra engine-specific config")
    args = ap.parse_args()

    _safe_mkdir(args.out_dir)

    df = pd.read_csv(args.input_csv)
    id_col = args.id_col or _pick_id_col(df, ["ID", "Name", "smiles", "MolID"])

    _require_col(df, args.top_col)
    _require_col(df, args.traj_col)

    cfg: Dict[str, Any] = {
        "stride": args.stride,
        "ref_frame": args.ref_frame,
        "sel_backbone": args.sel_backbone,
        "sel_protein": args.sel_protein,
        "sel_ligand": args.sel_ligand,
        "hb_donors_sel": args.hb_donors_sel,
        "hb_acceptors_sel": args.hb_acceptors_sel,
        "hb_distance": args.hb_distance,
        "hb_angle": args.hb_angle,
    }

    if args.config_json:
        with open(args.config_json, "r", encoding="utf-8") as f:
            extra = json.load(f)
        if isinstance(extra, dict):
            cfg.update(extra)

    fastmd = try_import_fastmdanalysis() if args.engine in ("auto", "fastmd") else None
    if args.engine == "fastmd" and fastmd is None:
        print("ERROR: engine=fastmd but fastmdanalysis cannot be imported.", file=sys.stderr)
        sys.exit(2)

    summaries: List[Dict[str, Any]] = []
    n_total = len(df)
    n_ok = 0
    n_skip = 0

    for i, row in df.iterrows():
        sys_id = str(row[id_col])
        top_path = str(row.get(args.top_col, ""))
        traj_path = str(row.get(args.traj_col, ""))

        if not (_file_exists(top_path) and _file_exists(traj_path)):
            summaries.append({
                "ID": sys_id,
                "status": "skipped",
                "reason": "missing top/traj file",
                "top_path": top_path,
                "traj_path": traj_path,
            })
            n_skip += 1
            print(f"[{i+1}/{n_total}] {sys_id}: skipped (missing files)")
            continue

        sys_out = os.path.join(args.out_dir, sys_id.replace("/", "_"))
        _safe_mkdir(sys_out)

        # Prefer fastmdanalysis if available; fallback to MDAnalysis if auto
        if fastmd is not None and args.engine in ("auto", "fastmd"):
            res = analyze_with_fastmd(fastmd, top_path, traj_path, sys_out, cfg)
            if res.get("status") != "ok" and args.engine == "auto":
                res["fallback_reason"] = "fastmd_failed"
                res = analyze_with_mdanalysis(top_path, traj_path, sys_out, cfg)
        else:
            res = analyze_with_mdanalysis(top_path, traj_path, sys_out, cfg)

        res["ID"] = sys_id
        summaries.append(res)

        if res.get("status") == "ok":
            n_ok += 1

        print(f"[{i+1}/{n_total}] {sys_id}: {res.get('engine')} -> {res.get('status')}")

    out_summary = os.path.join(args.out_dir, "step6a_md_summary.csv")
    pd.DataFrame(summaries).to_csv(out_summary, index=False)

    print(f"\n✅ Done. OK={n_ok}, skipped={n_skip}, total={n_total}")
    print(f"✅ Summary saved: {out_summary}")


if __name__ == "__main__":
    main()
