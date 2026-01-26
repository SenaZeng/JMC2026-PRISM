# -*- coding: utf-8 -*-
"""
core/step6c_run_gromacs.py  (Step6C - PREPARE ONLY)
==================================================
目标：为每个候选分子/replica 生成 MD 目录与“参数化准备脚本”run_md.sh（不直接跑 MD 生产段）。

本 Step6C 的产物（每个 job）应当稳定生成：
  md/<mol>/r<k>/01_param/complex.gro
  md/<mol>/r<k>/01_param/topol.top     (只包含 Protein + ligand 的唯一 [molecules] 段；不含 SOL/NA/CL)
  md/<mol>/r<k>/01_param/ligand.itp
  md/<mol>/r<k>/01_param/ligand.gro
  md/<mol>/r<k>/01_param/posre.itp
  md/<mol>/r<k>/01_param/ligand_moltype.txt

关键修复点（针对你遇到的失败模式）：
1) ACPYPE 参数化增加 CHARGE sweep（避免 odd electrons / sqm 失败导致 10/15 jobs 失败）；
   - 默认尝试：CHARGE(用户给定) + [0, 1, -1, 2, -2, 3, -3]
   - 可通过环境变量 CHARGE_LIST 覆盖，例如：CHARGE_LIST="0 1 -1"
2) slurm_array.sbatch 在 run_md.sh 失败时回显 ligand_acpype/acpype_charge_*.log 的尾部，便于定位。
3) 解析 ligand.itp 的 [ moleculetype ] 名称并写入 topol.top 的 [ molecules ]。
4) 强制删除 topol.top 中任何已有 [ molecules ] 段，仅保留一个正确段，避免重复导致拓扑错配。
5) ligand 提取：只用 HETATM（不 fallback ATOM），避免把蛋白 ATOM 当 ligand。
6) ACPYPE 输出路径：用 find 自动拾取 *_GMX.itp/gro，避免 glob 失败。
7) pdb2gmx 使用 protein_only.pdb + -ignh，避免输入含氢/UNL 出错。

用法：
  python core/step6c_run_gromacs.py --prepare_only --scheduler slurm --slurm_partition cpuPartition
  sbatch md/slurm_array.sbatch
"""

import os
import argparse
from pathlib import Path
import pandas as pd


# ----------------- Path setup ----------------- #
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))

DEFAULT_MD_INPUT = os.path.join(PROJECT_ROOT, "results", "step6b_md_input.csv")
DEFAULT_MD_DIR = os.path.join(PROJECT_ROOT, "md")


def write_text(p: Path, text: str) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


RUN_MD_SH_TEMPLATE = r"""#!/usr/bin/env bash
set -euo pipefail

# -------------- User config --------------
GMX=${GMX:-gmx}
FF_PROTEIN=${FF_PROTEIN:-amber99sb-ildn}
WATER_MODEL=${WATER_MODEL:-tip3p}

USE_ACPYPE=${USE_ACPYPE:-1}      # 1: run acpype from ligand.pdb
CHARGE=${CHARGE:-0}              # preferred net charge for acpype -n (will be tried first)
# Optional override: export CHARGE_LIST="0 1 -1 2 -2" (space-separated)
CHARGE_LIST=${CHARGE_LIST:-}

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-__OMP_THREADS__}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IN_DIR="${ROOT_DIR}/00_inputs"
PARAM_DIR="${ROOT_DIR}/01_param"

mkdir -p "$PARAM_DIR"
cd "$PARAM_DIR"
cp "${IN_DIR}/complex.pdb" .

echo "[PREP] ROOT_DIR=${ROOT_DIR}"
echo "[PREP] GMX=$($GMX --version | head -n 1)"
echo "[PREP] OMP_NUM_THREADS=${OMP_NUM_THREADS}"

# ---------------- 1) protein_only.pdb + pdb2gmx ----------------
python - << 'PY'
from pathlib import Path
lines = Path("complex.pdb").read_text().splitlines()
out = [l for l in lines if l.startswith("ATOM") or l.startswith("TER")]
out.append("END")
Path("protein_only.pdb").write_text("\n".join(out) + "\n")
print(f"[protein_only] atoms: {sum(1 for l in out if l.startswith('ATOM'))}")
PY

# Choose default termini (group 1) non-interactively
$GMX pdb2gmx -f protein_only.pdb -ignh \
  -o protein.gro -p topol.top -i posre.itp \
  -ff "$FF_PROTEIN" -water "$WATER_MODEL" << EOF
1
EOF

# ---------------- 2) Extract ligand (HETATM only) ----------------
python - << 'PY'
from pathlib import Path
inp = Path("complex.pdb").read_text().splitlines()
lig = [l for l in inp if l.startswith("HETATM")]
if not lig:
    raise SystemExit("[ERR] No HETATM found for ligand in complex.pdb (cannot continue).")
Path("ligand.pdb").write_text("\n".join(lig) + "\nEND\n")
print(f"[ligand] atoms: {len(lig)}")
PY

# ---------------- 3) Ligand parameterization (ACPYPE/GAFF) ----------------
if [[ "$USE_ACPYPE" == "1" ]]; then
  rm -rf ligand_acpype
  mkdir -p ligand_acpype
  cp ligand.pdb ligand_acpype/

  # Build charge candidates:
  # - If CHARGE_LIST provided, use it as-is.
  # - Else try: user CHARGE first, then common charges.
  if [[ -n "$CHARGE_LIST" ]]; then
    read -r -a _CHARGE_CANDIDATES <<< "$CHARGE_LIST"
  else
    _CHARGE_CANDIDATES=("$CHARGE" 0 1 -1 2 -2 3 -3)
  fi

  # de-dup safely (bash 4.2 compatible)
  _ACPYPE_OK=0
  _ACPYPE_USED_CHARGE=""

  for q in "${_CHARGE_CANDIDATES[@]}"; do
    echo "[ACPYPE] try CHARGE=$q"
    rm -f ligand_acpype/acpype_charge_*.log 2>/dev/null || true

    # Keep per-charge log
    (cd ligand_acpype && acpype -i ligand.pdb -b ligand -n "$q" > "acpype_charge_${q}.log" 2>&1) || true

    # robust pickup of ACPYPE outputs
    LIG_ITP=$(find ligand_acpype -maxdepth 4 -name "*_GMX.itp" | head -n 1 || true)
    LIG_GRO=$(find ligand_acpype -maxdepth 4 -name "*_GMX.gro" | head -n 1 || true)

    if [[ -n "$LIG_ITP" && -n "$LIG_GRO" ]]; then
      _ACPYPE_OK=1
      _ACPYPE_USED_CHARGE="$q"
      break
    fi

    echo "[ACPYPE] CHARGE=$q failed (no *_GMX.itp/gro)."
  done

  if [[ "$_ACPYPE_OK" != "1" ]]; then
    echo "[ERR] ACPYPE failed for all charges tried: ${_CHARGE_CANDIDATES[*]}"
    echo "[ERR] Logs saved under ligand_acpype/acpype_charge_*.log"
    echo "[ERR] Tree:"
    find ligand_acpype -maxdepth 4 -type f | head -n 200
    exit 3
  fi



  cp "$LIG_ITP" ligand.itp
  cp "$LIG_GRO" ligand.gro
  echo "[ACPYPE] success with CHARGE=$_ACPYPE_USED_CHARGE"
  echo "[ACPYPE] picked itp=$LIG_ITP"
  echo "[ACPYPE] picked gro=$LIG_GRO"
else
  echo "[ERR] USE_ACPYPE=0 is not supported in this pipeline yet."
  exit 4
fi

# ---------------- 4) Merge protein + ligand coordinates -> complex.gro ----------------
$GMX editconf -f protein.gro -o protein.pdb

python - << 'PY'
from pathlib import Path

prot = Path("protein.pdb").read_text().splitlines()
gro = Path("ligand.gro").read_text().splitlines()

nat = int(gro[1].strip())
atoms = gro[2:2+nat]  # fixed-format .gro atom lines

out = [l for l in prot if not l.startswith("END")]

# Convert nm -> Å, write placeholder PDB HETATM lines (enough for editconf)
i = 1
for a in atoms:
    x = float(a[20:28]) * 10.0
    y = float(a[28:36]) * 10.0
    z = float(a[36:44]) * 10.0
    out.append(f"HETATM{i:5d}  C   LIG A9999    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
    i += 1

out.append("END")
Path("complex_merged.pdb").write_text("\n".join(out) + "\n")
PY

$GMX editconf -f complex_merged.pdb -o complex.gro

# ---------------- 5) Make topol.top include ligand + SINGLE [molecules] (NO SOL) ----------------
python - << 'PY'
from pathlib import Path
import re

def parse_lig_moltype(itp_path: Path) -> str:
    lines = itp_path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.strip().lower() == "[ moleculetype ]":
            for j in range(i+1, min(i+50, len(lines))):
                s = lines[j].strip()
                if (not s) or s.startswith(";") or s.startswith("#"):
                    continue
                return s.split()[0]
    raise RuntimeError("Cannot parse ligand moleculetype from ligand.itp")

top_path = Path("topol.top")
itp_path = Path("ligand.itp")

lig_name = parse_lig_moltype(itp_path)
Path("ligand_moltype.txt").write_text(lig_name + "\n")
print(f"[TOPO] ligand moleculetype = {lig_name}")

# Read current topol.top
out = top_path.read_text().splitlines()

include_line = '#include "ligand.itp"'

# 1) remove any existing ligand include lines (dedupe; robust to spaces/quotes)
out = [ln for ln in out if not re.match(r'^\s*#include\s+["\']ligand\.itp["\']\s*$', ln)]

# 2) insert include right after forcefield include if possible; otherwise insert at top
idx_ff = None
for i, ln in enumerate(out):
    if "forcefield.itp" in ln:
        idx_ff = i
        break

if idx_ff is not None:
    out.insert(idx_ff + 1, include_line)
else:
    out.insert(0, include_line)

# 3) remove ALL existing [ molecules ] blocks (prevents duplicates)
clean = []
i = 0
while i < len(out):
    if re.match(r'^\s*\[\s*molecules\s*\]\s*$', out[i], re.I):
        # skip until next directive or EOF
        i += 1
        while i < len(out) and not re.match(r'^\s*\[.*\]\s*$', out[i]):
            i += 1
        continue
    clean.append(out[i])
    i += 1

out = clean

# 3) append exactly ONE molecules block, and DO NOT include SOL/ions here
clean.append("")
clean.append("[ molecules ]")
clean.append("; Compound        #mols")
clean.append("Protein_chain_A    1")
clean.append(f"{lig_name:<16s} 1")
clean.append("")

top_path.write_text("\n".join(clean) + "\n")
print("[TOPO] topol.top rewritten with single [molecules] (Protein + ligand only).")
PY

echo "[OK] PREP DONE: generated 01_param/{complex.gro,topol.top,ligand.itp,ligand.gro,posre.itp}"
"""

SLURM_ARRAY_TEMPLATE = r"""#!/usr/bin/env bash
#SBATCH -J mpro_prep
#SBATCH -p __PARTITION__
#SBATCH --cpus-per-task=__CPUS__
#SBATCH --mem=__MEM__
#SBATCH -t __TIME__
#SBATCH --array=0-__ARRAY_MAX__

set -euo pipefail

# -----------------------------
# FIX: conda activate/deactivate scripts may reference undefined vars,
# which will crash under "set -u". Also avoid inheriting login-shell conda state.
# -----------------------------
unset CONDA_SHLVL CONDA_PREFIX CONDA_DEFAULT_ENV

# Temporarily disable nounset during conda init/activate
set +u
source ~/anaconda3/etc/profile.d/conda.sh
conda activate __CONDA_ENV__
set -u

export OMP_NUM_THREADS=__CPUS__

JOBS=(
__JOBS_LIST__
)

JOB="${JOBS[$SLURM_ARRAY_TASK_ID]}"
JOB_DIR="$(dirname "$JOB")"
echo "[SLURM] Running: $JOB"

set +e
(cd "$JOB_DIR" && bash "$(basename "$JOB")")
RET=$?
set -e

if [[ $RET -ne 0 ]]; then
  echo "================== [SLURM][ERROR] run_md.sh failed (exit=$RET) =================="
  echo "[SLURM][ERROR] JOB=$JOB"
  echo "[SLURM][ERROR] PWD=$(pwd)"
  if [[ -d "$JOB_DIR/ligand_acpype" ]]; then
    echo "---- ligand_acpype logs (tail -n 120) ----"
    ls -lah "$JOB_DIR/ligand_acpype" || true
    tail -n 120 "$JOB_DIR/ligand_acpype"/acpype_charge_*.log 2>/dev/null || true
  else
    echo "[SLURM][ERROR] $JOB_DIR/ligand_acpype not found"
  fi
  echo "==========================================================================="
  exit $RET
fi
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--md_input", type=str, default=DEFAULT_MD_INPUT)
    ap.add_argument("--out_md_dir", type=str, default=DEFAULT_MD_DIR)

    ap.add_argument("--prepare_only", action="store_true", help="Only generate scripts/dirs; do not run.")
    ap.add_argument("--scheduler", type=str, default="none", choices=["none", "slurm"])

    # Slurm options
    ap.add_argument("--slurm_partition", type=str, default="cpu")
    ap.add_argument("--slurm_cpus", type=int, default=8)
    ap.add_argument("--slurm_mem", type=str, default="16G")
    ap.add_argument("--slurm_time", type=str, default="48:00:00")
    ap.add_argument("--conda_env", type=str, default="md_env")

    args = ap.parse_args()

    md_input = args.md_input if os.path.isabs(args.md_input) else os.path.join(PROJECT_ROOT, args.md_input)
    out_md_dir = args.out_md_dir if os.path.isabs(args.out_md_dir) else os.path.join(PROJECT_ROOT, args.out_md_dir)

    print(f"[INFO] project_root = {PROJECT_ROOT}")
    print(f"[INFO] md_input     = {md_input}")
    print(f"[INFO] out_md_dir   = {out_md_dir}")

    df = pd.read_csv(md_input)
    out_md_path = Path(out_md_dir)
    out_md_path.mkdir(parents=True, exist_ok=True)

    jobs = []

    # Build md/<mol>/r<k>/00_inputs/complex.pdb and run_md.sh
    for _, row in df.iterrows():
        name = str(row["name"])
        replicas = int(row.get("replicas", 1))

        complex_pdb = str(row["complex_pdb"])
        if not os.path.isabs(complex_pdb):
            complex_pdb = os.path.join(PROJECT_ROOT, complex_pdb)

        for k in range(1, replicas + 1):
            job_root = out_md_path / name / f"r{k}"
            in_dir = job_root / "00_inputs"
            in_dir.mkdir(parents=True, exist_ok=True)

            # Copy complex.pdb into 00_inputs
            dst_complex = in_dir / "complex.pdb"
            if not dst_complex.exists():
                dst_complex.write_bytes(Path(complex_pdb).read_bytes())

            # Write run_md.sh (prepare-only pipeline)
            run_sh = job_root / "run_md.sh"
            omp = args.slurm_cpus if args.scheduler == "slurm" else 8
            script = RUN_MD_SH_TEMPLATE.replace("__OMP_THREADS__", str(omp))
            write_text(run_sh, script)
            os.chmod(run_sh, 0o755)

            jobs.append(str(run_sh))

    # local runner
    local_runner = out_md_path / "run_all_local.sh"
    local_txt = "#!/usr/bin/env bash\nset -euo pipefail\n\n"
    local_txt += "JOBS=(\n" + "\n".join([f'  "{p}"' for p in jobs]) + "\n)\n\n"
    local_txt += 'for f in "${JOBS[@]}"; do (cd "$(dirname "$f")" && bash "$(basename "$f")"); done\n'
    write_text(local_runner, local_txt)
    os.chmod(local_runner, 0o755)

    # Slurm array script
    if args.scheduler == "slurm":
        slurm = out_md_path / "slurm_array.sbatch"
        jobs_list = "\n".join([f'  "{p}"' for p in jobs])
        sb = SLURM_ARRAY_TEMPLATE
        sb = sb.replace("__PARTITION__", args.slurm_partition)
        sb = sb.replace("__CPUS__", str(args.slurm_cpus))
        sb = sb.replace("__MEM__", args.slurm_mem)
        sb = sb.replace("__TIME__", args.slurm_time)
        sb = sb.replace("__CONDA_ENV__", args.conda_env)
        sb = sb.replace("__ARRAY_MAX__", str(max(len(jobs) - 1, 0)))
        sb = sb.replace("__JOBS_LIST__", jobs_list)
        write_text(slurm, sb)

    print(f"\n[DONE] Prepared jobs = {len(jobs)}")
    print(f"[DONE] Local runner: {local_runner}")
    if args.scheduler == "slurm":
        print(f"[DONE] Slurm script: {out_md_path / 'slurm_array.sbatch'}")
    print("\nNext:")
    print("  sbatch md/slurm_array.sbatch    # run PREP for all jobs (generate 01_param outputs)")


if __name__ == "__main__":
    main()