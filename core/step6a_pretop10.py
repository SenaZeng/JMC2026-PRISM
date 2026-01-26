# -*- coding: utf-8 -*-
"""
Step6A - MD Pre-selection (Top 5-10)  [Accumulate + Full Trace Columns + Report]
--------------------------------------------------------------------------------
目标：为 Docking → MD 的“工程断层”提供节流阀，并保证 TopN 尽量不为空/不全灭。

筛选机制（核心思想）：
A) Hard Filter（物理/任务约束）：
   1) 可选：Filter_Status == Pass（--require_pass 开启时）
   2) Broad_Spectrum_Score <= broad_thr（越负越好，代表“最差靶点也强”的广谱能力）
   3) PySCF_Gap_eV 落入给定 gap_windows 的窗口集合（先黄金窗，再逐步扩窗补齐）

B) Accumulate to TopN（避免全灭/不足）：
   - 按 gap_windows 从窄到宽逐个筛选；
   - 每个窗口产出的候选“追加”到 selected；
   - 按 dedup_key（name 或 smiles）去重；
   - 每次追加后按统一 ranking 重新排序；
   - 直到 selected >= TopN 或窗口用尽。

C) Ranking（只决定顺序，不决定是否入池）：
   主：Broad_Spectrum_Score（越负越好）
   次：gap_dev = |PySCF_Gap_eV - 5.0|（越接近 5 eV 越好）
   再：R_total（若存在，越大越好，作为 tie-break）

输出：
- results/step6a_md_gold_candidates_top.csv（TopN 榜单 + 追溯列）
- results/step6a_md_joblist.csv（MD 批量任务单）
- results/step6a_md_prescreen_report.txt（筛选机制说明，可直接贴报告）

默认路径（project_root/results/...）：
- in final:  results/step5b_final_candidates.csv
- out top:   results/step6a_md_gold_candidates_top.csv
- out jobs:  results/step6a_md_joblist.csv
"""

import os
import argparse
from typing import Tuple, Optional, List

import pandas as pd


# ----------------- 路径基础设置 ----------------- #
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))

DEFAULT_IN_FINAL_5B = os.path.join(project_root, "results", "step5b_final_candidates.csv")
DEFAULT_OUT_TOP = os.path.join(project_root, "results", "step6a_md_gold_candidates_top.csv")
DEFAULT_OUT_JOBLIST = os.path.join(project_root, "results", "step6a_md_joblist.csv")
DEFAULT_OUT_REPORT = os.path.join(project_root, "results", "step6a_md_prescreen_report.txt")


def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    if "Filter_Status" not in out.columns:
        out["Filter_Status"] = "Pass"
    out["Filter_Status"] = out["Filter_Status"].astype(str)

    if "R_total" not in out.columns and "R_global" in out.columns:
        out["R_total"] = out["R_global"]

    required = ["Broad_Spectrum_Score", "PySCF_Gap_eV"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    if "name" not in out.columns:
        out["name"] = ""
    if "smiles" not in out.columns:
        out["smiles"] = ""

    return out


def _rank_sort(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    sort_cols = ["Broad_Spectrum_Score", "gap_dev"]
    ascending = [True, True]

    if "R_total" in df.columns:
        sort_cols.append("R_total")
        ascending.append(False)

    return df.sort_values(by=sort_cols, ascending=ascending).reset_index(drop=True)


def _prescreen_pool(
    df: pd.DataFrame,
    broad_thr: float,
    gap_window: Tuple[float, float],
    require_pass: bool = True,
) -> pd.DataFrame:
    glo, ghi = gap_window

    cond = (df["Broad_Spectrum_Score"] <= broad_thr) & df["PySCF_Gap_eV"].between(glo, ghi)
    if require_pass:
        cond = cond & (df["Filter_Status"].str.lower() == "pass")

    pool = df.loc[cond].copy()
    if pool.empty:
        return pool

    pool["gap_dev"] = (pool["PySCF_Gap_eV"] - 5.0).abs()
    return _rank_sort(pool)


def _accumulate_to_topn(
    df: pd.DataFrame,
    broad_thr: float,
    gap_windows: List[Tuple[float, float]],
    top_n: int,
    require_pass: bool,
    dedup_key: str,
) -> Tuple[pd.DataFrame, str]:
    selected = pd.DataFrame()
    used_windows: List[Tuple[float, float]] = []

    for w in gap_windows:
        pool_w = _prescreen_pool(df=df, broad_thr=broad_thr, gap_window=w, require_pass=require_pass)
        if pool_w.empty:
            continue

        used_windows.append(w)

        if dedup_key not in pool_w.columns:
            raise ValueError(f"dedup_key '{dedup_key}' not found in columns")

        if not selected.empty:
            pool_w = pool_w[~pool_w[dedup_key].isin(selected[dedup_key])].copy()

        if pool_w.empty:
            continue

        selected = pd.concat([selected, pool_w], ignore_index=True)
        selected = _rank_sort(selected)

        if len(selected) >= top_n:
            break

    used_window_str = ";".join([f"{lo}-{hi}" for lo, hi in used_windows]) if used_windows else ""
    return selected, used_window_str


def _build_joblist(top: pd.DataFrame, engine: str, md_length_ns: int, replicas_top3: int) -> pd.DataFrame:
    job = top[["rank", "name", "smiles"]].copy()
    job["engine"] = engine
    job["md_length_ns"] = int(md_length_ns)
    job["replicas"] = 1
    job.loc[job["rank"] <= 3, "replicas"] = int(replicas_top3)
    return job


def _write_report(path: str, args, used_windows: str, n_selected: int) -> None:
    txt = f"""[Step6A] MD 预选机制说明（可直接贴报告）

输入：
- Step5B final candidates: {args.in_final}

硬过滤（Hard Filter）：
1) require_pass_used = {bool(args.require_pass)}
   - 若开启：仅保留 Filter_Status == "Pass" 的分子进入候选池；
2) broad_thr_used = {args.broad_thr}
   - 仅保留 Broad_Spectrum_Score <= broad_thr 的分子（越负越好，强调“最差靶点也强”的广谱短板能力）；
3) gap_windows_used = {used_windows if used_windows else "FALLBACK(no_gap_filter)"}
   - 先用较窄“黄金窗”筛选，若数量不足 TopN，则逐步扩窗补齐；
4) dedup_key_used = {args.dedup_key}
   - 在不同 gap 窗口累计补位时，按该字段去重（防止同一分子重复进入候选池）。

排序（Ranking，仅影响顺序）：
- 主排序：Broad_Spectrum_Score（越负越好）
- 次排序：gap_dev = |PySCF_Gap_eV - 5.0|（越接近 5 eV 越好）
- 追加排序：R_total（若存在，越大越好，作为 tie-break）

输出：
- TopN 榜单：{args.out_top}
- MD 任务单：{args.out_joblist}

备注：
- 本步骤用于“节流”：优先把算力投入到广谱强且电子结构更稳定（Gap 更接近理想值）的分子；
- 若强约束导致数量不足，将通过 gap_windows 的扩窗策略补齐（仍保持 broad_thr 与 require_pass 的核心标准）。
"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(txt)


def main():
    parser = argparse.ArgumentParser(description="Step6A: MD Pre-selection TopK (Accumulate + Trace)")
    parser.add_argument("--in_final", type=str, default=DEFAULT_IN_FINAL_5B)
    parser.add_argument("--out_top", type=str, default=DEFAULT_OUT_TOP)
    parser.add_argument("--out_joblist", type=str, default=DEFAULT_OUT_JOBLIST)
    parser.add_argument("--out_report", type=str, default=DEFAULT_OUT_REPORT)

    parser.add_argument("--top_n", type=int, default=10)
    parser.add_argument("--broad_thr", type=float, default=-7.4)
    parser.add_argument(
        "--gap_windows",
        type=str,
        default="4.7,5.3;4.6,5.4;4.5,5.5;4.4,5.6",
        help="format: 'lo,hi;lo,hi;...' (ordered, first is preferred)",
    )
    parser.add_argument("--require_pass", action="store_true")
    parser.add_argument("--dedup_key", type=str, default="name", choices=["name", "smiles"])

    parser.add_argument("--engine", type=str, default="gromacs", choices=["gromacs", "amber"])
    parser.add_argument("--md_length_ns", type=int, default=100)
    parser.add_argument("--replicas_top3", type=int, default=3)

    args = parser.parse_args()

    if not os.path.exists(args.in_final):
        raise FileNotFoundError(f"Input not found: {args.in_final}")

    df = pd.read_csv(args.in_final)
    df = _normalize_columns(df)

    # parse gap windows
    gap_windows: List[Tuple[float, float]] = []
    for part in args.gap_windows.split(";"):
        lo_str, hi_str = part.split(",")
        gap_windows.append((float(lo_str.strip()), float(hi_str.strip())))
    if not gap_windows:
        raise ValueError("No valid gap_windows parsed")

    pool, used_window_str = _accumulate_to_topn(
        df=df,
        broad_thr=args.broad_thr,
        gap_windows=gap_windows,
        top_n=args.top_n,
        require_pass=bool(args.require_pass),
        dedup_key=args.dedup_key,
    )

    # fallback: guarantee non-empty output
    if pool.empty:
        fallback = df.copy()
        fallback["gap_dev"] = (fallback["PySCF_Gap_eV"] - 5.0).abs()
        pool = _rank_sort(fallback)
        used_window_str = "FALLBACK(no_gap_filter)"

    top = pool.head(args.top_n).copy()
    top.insert(0, "rank", range(1, len(top) + 1))

    # trace columns (always present in output)
    top["gap_windows_used"] = used_window_str
    top["broad_thr_used"] = args.broad_thr
    top["require_pass_used"] = bool(args.require_pass)
    top["dedup_key_used"] = args.dedup_key

    # output top
    os.makedirs(os.path.dirname(args.out_top), exist_ok=True)

    cols = [
        "rank", "name", "smiles",
        "Broad_Spectrum_Score", "PySCF_Gap_eV", "gap_dev"
    ]
    if "R_total" in top.columns:
        cols.append("R_total")

    # 强制展示你想要的追溯列
    cols += ["gap_windows_used", "broad_thr_used", "require_pass_used", "dedup_key_used"]

    top[cols].to_csv(args.out_top, index=False)

    # output joblist
    os.makedirs(os.path.dirname(args.out_joblist), exist_ok=True)
    job = _build_joblist(top=top, engine=args.engine, md_length_ns=args.md_length_ns, replicas_top3=args.replicas_top3)
    job.to_csv(args.out_joblist, index=False)

    # report
    _write_report(args.out_report, args, used_window_str, len(top))

    # print: avoid pandas truncation
    pd.set_option("display.max_columns", 200)
    pd.set_option("display.width", 200)

    print(f"[Step6A] Input:  {args.in_final}")
    print(f"[Step6A] Output: {args.out_top}")
    print(f"[Step6A] Jobs:   {args.out_joblist}")
    print(f"[Step6A] Report: {args.out_report}")
    print(top[cols].to_string(index=False))


if __name__ == "__main__":
    main()
