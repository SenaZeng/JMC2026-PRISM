#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step6d_run_gromacs_run.py (replacement)
- Runs GROMACS MD for each (mol_x, replica) prepared by step6c/step6b.
- Adds automatic QC + fallback ("rescue") between EM(min) and NVT.

Workflow per job:
01_param (from step6c)  -> 03_min (EM) -> QC -> 04_nvt (main or fallback) -> 05_npt -> 06_md

QC metrics (after EM):
- Fmax parsed from em.log ("Maximum force")
- Minimum distance (System vs SOL) from gmx mindist on em.gro

If QC fails:
- strong EM (emtol lower + more steps)
- short 50K NVT (20 ps)
- 300K NVT (80-100 ps)
then continue with NPT -> MD
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


@dataclass
class Job:
    name: str
    replica: int
    md_ns: float
    root: Path  # md/{name}/r{replica}
    param_dir: Path  # root/01_param


def _read_jobs(md_input_csv: Path, md_root: Path) -> List[Job]:
    jobs: List[Job] = []
    with md_input_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"].strip()
            reps = int(float(row.get("replicas", 1) or 1))
            md_ns = float(row.get("md_length_ns", 100) or 100)
            for k in range(1, reps + 1):
                root = md_root / name / f"r{k}"
                param_dir = root / "01_param"
                jobs.append(Job(name=name, replica=k, md_ns=md_ns, root=root, param_dir=param_dir))
    return jobs


def _sh_escape(s: str) -> str:
    return "'" + s.replace("'", "'\"'\"'") + "'"


def render_run_md2_sh(
    job: Job,
    threads: int,
    mindist_thr_nm: float,
    fmax_thr: float,
    fb_emtol: float,
    fb_em_nsteps: int,
    fb_nvt50_ps: float,
    fb_nvt300_ps: float,
    main_nvt_ps: float,
    npt_ps: float,
) -> str:
    root = str(job.root.resolve())
    md_ps = job.md_ns * 1000.0
    # NOTE: We always apply POSRES during NVT/NPT equilibration stages.
    # Production MD has no POSRES.
    return f"""#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR={_sh_escape(root)}
GMX="${{GMX:-gmx}}"
export OMP_NUM_THREADS="${{OMP_NUM_THREADS:-{threads}}}"

MD_NS={job.md_ns}
MAIN_NVT_PS={main_nvt_ps}
NPT_PS={npt_ps}

QC_MINDIST_THR_NM={mindist_thr_nm}
QC_FMAX_THR={fmax_thr}

FB_EMTOL={fb_emtol}
FB_EM_NSTEPS={fb_em_nsteps}
FB_NVT50_PS={fb_nvt50_ps}
FB_NVT300_PS={fb_nvt300_ps}

echo "[MD2] ROOT_DIR=$ROOT_DIR"
echo "[MD2] GMX=$($GMX --version | head -n 1 || true)"
echo "[MD2] OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "[MD2] MD_NS=$MD_NS  MAIN_NVT_PS=$MAIN_NVT_PS  NPT_PS=$NPT_PS"
echo "[MD2] QC: mindist_thr=$QC_MINDIST_THR_NM nm  fmax_thr=$QC_FMAX_THR"
echo "[MD2] FB: emtol=$FB_EMTOL nsteps=$FB_EM_NSTEPS nvt50=${{FB_NVT50_PS}}ps nvt300=${{FB_NVT300_PS}}ps"

cd "$ROOT_DIR"

PARAM_DIR="01_param"
if [[ ! -d "$PARAM_DIR" ]]; then
  echo "[MD2][ERR] Missing $PARAM_DIR (step6c not prepared?)"
  exit 2
fi

# Convenience copies
TOPOL="$PARAM_DIR/topol.top"
LIGITP="$PARAM_DIR/ligand.itp"
POSRE="$PARAM_DIR/posre.itp"

# ---------- helpers ----------
extract_fmax_from_emlog() {{
  local emlog="$1"
  # line example: "Maximum force     =  8.0680969e+01 on atom 60881"
  awk '
    $1=="Maximum" && $2=="force" {{
      for(i=1;i<=NF;i++) if($i=="=") print $(i+1)
    }}
  ' "$emlog" | tail -n 1
}}

calc_mindist_nm() {{
  local gro="$1"
  local tpr="$2"
  local outxvg="$3"
  # Choose groups non-interactively: 0(System) then 16(SOL) in your index ordering.
  # If your SOL group number differs, mindist will fail; then we fall back to System vs Water_and_ions (19) if available.
  if echo -e "0\\n16\\n" | $GMX mindist -f "$gro" -s "$tpr" -od "$outxvg" >/dev/null 2>&1; then
    tail -n 1 "$outxvg" | awk '{{print $2}}'
    return 0
  fi
  if echo -e "0\\n19\\n" | $GMX mindist -f "$gro" -s "$tpr" -od "$outxvg" >/dev/null 2>&1; then
    tail -n 1 "$outxvg" | awk '{{print $2}}'
    return 0
  fi
  echo "nan"
  return 0
}}

needs_repair() {{
  local fmax="$1"
  local mindist="$2"
  python - "$fmax" "$mindist" "$QC_FMAX_THR" "$QC_MINDIST_THR_NM" <<'PY'
import sys, math

fmax_s, mind_s, thr_f_s, thr_m_s = sys.argv[1:5]

def to_float(x):
    x = str(x).strip()
    if x in ("", "nan", "NaN", "None"):
        return float("nan")
    return float(x)

fmax = to_float(fmax_s)
mind = to_float(mind_s)
thr_f = to_float(thr_f_s)
thr_m = to_float(thr_m_s)

bad = False
if not math.isnan(mind) and not math.isnan(thr_m) and mind < thr_m:
    bad = True
if not math.isnan(fmax) and not math.isnan(thr_f) and fmax > thr_f:
    bad = True

print("1" if bad else "0")
PY
}}

write_mdp_em() {{
  local mdp="$1"
  local emtol="$2"
  local nsteps="$3"
  cat > "$mdp" <<EOF
integrator  = steep
emtol       = $emtol
emstep      = 0.01
nsteps      = $nsteps
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
EOF
}}

write_mdp_nvt() {{
  local mdp="$1"
  local nsteps="$2"
  local dt="$3"
  local ref_t="$4"
  cat > "$mdp" <<EOF
define = -DPOSRES
integrator  = md
nsteps      = $nsteps
dt          = $dt
continuation = no
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.2
ref_t       = $ref_t
pcoupl      = no
pbc         = xyz
gen_vel     = yes
gen_temp    = $ref_t
gen_seed    = -1
nstxout     = 0
nstvout     = 0
nstenergy   = 500
nstlog      = 500
EOF
}}

write_mdp_nvt_continue() {{
  local mdp="$1"
  local nsteps="$2"
  local dt="$3"
  local ref_t="$4"
  cat > "$mdp" <<EOF
define = -DPOSRES
integrator  = md
nsteps      = $nsteps
dt          = $dt
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.2
ref_t       = $ref_t
pcoupl      = no
pbc         = xyz
gen_vel     = no
nstxout     = 0
nstvout     = 0
nstenergy   = 500
nstlog      = 500
EOF
}}

write_mdp_npt_continue() {{
  local mdp="$1"
  local nsteps="$2"
  local dt="$3"
  local ref_t="$4"
  cat > "$mdp" <<EOF
define = -DPOSRES
integrator  = md
nsteps      = $nsteps
dt          = $dt
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.2
ref_t       = $ref_t
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
gen_vel     = no
nstxout     = 0
nstvout     = 0
nstenergy   = 500
nstlog      = 500
EOF
}}

write_mdp_md() {{
  local mdp="$1"
  local nsteps="$2"
  local dt="$3"
  local ref_t="$4"
  cat > "$mdp" <<EOF
integrator  = md
nsteps      = $nsteps
dt          = $dt
continuation = yes
constraint_algorithm = lincs
constraints = h-bonds
lincs_iter  = 2
lincs_order = 4
cutoff-scheme = Verlet
nstlist     = 20
rcoulomb    = 1.0
rvdw        = 1.0
coulombtype = PME
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.5
ref_t       = $ref_t
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
nstxout-compressed = 5000
compressed-x-grps   = System
nstxout     = 0
nstvout     = 0
nstenergy   = 1000
nstlog      = 1000
EOF
}}

# ---------- 03_min (EM) ----------
mkdir -p 03_min
cd 03_min

cp -f "../$PARAM_DIR/complex.gro" complex.gro
cp -f "../$PARAM_DIR/topol.top" topol.top
cp -f "../$PARAM_DIR/ligand.itp" ligand.itp
cp -f "../$PARAM_DIR/posre.itp" posre.itp

write_mdp_em em.mdp 500.0 50000
$GMX grompp -f em.mdp -c complex.gro -p topol.top -o em.tpr -maxwarn 5
$GMX mdrun -deffnm em -ntmpi 1 -ntomp "$OMP_NUM_THREADS"

FMAX="$(extract_fmax_from_emlog em.log || true)"
MINDIST="$({{ calc_mindist_nm em.gro em.tpr mindist_em.xvg; }} || true)"

echo "[MD2] QC after EM: FMAX=$FMAX  MINDIST_NM=$MINDIST"

DO_REPAIR="$(needs_repair "$FMAX" "$MINDIST")"
if [[ "$DO_REPAIR" == "1" ]]; then
  echo "[MD2][QC] needs_repair=YES -> running fallback (strong EM + staged NVT)"
  write_mdp_em em_strong.mdp "$FB_EMTOL" "$FB_EM_NSTEPS"
  $GMX grompp -f em_strong.mdp -c em.gro -p topol.top -o em_strong.tpr -maxwarn 5
  $GMX mdrun -deffnm em_strong -ntmpi 1 -ntomp "$OMP_NUM_THREADS"
  EM_GRO="em_strong.gro"
else
  echo "[MD2][QC] needs_repair=NO -> main path"
  EM_GRO="em.gro"
fi

cd ..

# ---------- 04_nvt ----------
mkdir -p 04_nvt
cd 04_nvt
cp -f "../03_min/$EM_GRO" start.gro
cp -f "../03_min/topol.top" topol.top
cp -f "../03_min/ligand.itp" ligand.itp
cp -f "../03_min/posre.itp" posre.itp

DT=0.002

if [[ "$DO_REPAIR" == "1" ]]; then
  # 50K, 20 ps (default)
  NSTEPS_50K=$(python - <<PY
print(int(round($FB_NVT50_PS / $DT)))
PY
)
  write_mdp_nvt nvt50.mdp "$NSTEPS_50K" "$DT" 50
  $GMX grompp -f nvt50.mdp -c start.gro -r start.gro -p topol.top -o nvt50.tpr -maxwarn 5
  $GMX mdrun -deffnm nvt50 -ntmpi 1 -ntomp "$OMP_NUM_THREADS"
  # 300K, 80 ps (default)
  NSTEPS_300K=$(python - <<PY
print(int(round($FB_NVT300_PS / $DT)))
PY
)
  write_mdp_nvt_continue nvt300.mdp "$NSTEPS_300K" "$DT" 300
  $GMX grompp -f nvt300.mdp -c nvt50.gro -r nvt50.gro -p topol.top -o nvt300.tpr -maxwarn 5
  $GMX mdrun -deffnm nvt300 -ntmpi 1 -ntomp "$OMP_NUM_THREADS"
  NVT_GRO="nvt300.gro"
else
  # main NVT at 300K
  NSTEPS_MAIN=$(python - <<PY
print(int(round($MAIN_NVT_PS / $DT)))
PY
)
  write_mdp_nvt nvt.mdp "$NSTEPS_MAIN" "$DT" 300
  $GMX grompp -f nvt.mdp -c start.gro -r start.gro -p topol.top -o nvt.tpr -maxwarn 5
  $GMX mdrun -deffnm nvt -ntmpi 1 -ntomp "$OMP_NUM_THREADS"
  NVT_GRO="nvt.gro"
fi

cd ..

# ---------- 05_npt ----------
mkdir -p 05_npt
cd 05_npt
cp -f "../04_nvt/$NVT_GRO" nvt.gro
cp -f "../04_nvt/topol.top" topol.top
cp -f "../04_nvt/ligand.itp" ligand.itp
cp -f "../04_nvt/posre.itp" posre.itp

NSTEPS_NPT=$(python - <<PY
print(int(round($NPT_PS / $DT)))
PY
)
write_mdp_npt_continue npt.mdp "$NSTEPS_NPT" "$DT" 300
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 5
$GMX mdrun -deffnm npt -ntmpi 1 -ntomp "$OMP_NUM_THREADS"

cd ..

# ---------- 06_md ----------
mkdir -p 06_md
cd 06_md
cp -f "../05_npt/npt.gro" npt.gro
cp -f "../05_npt/topol.top" topol.top
cp -f "../05_npt/ligand.itp" ligand.itp

NSTEPS_MD=$(python - <<PY
print(int(round(({md_ps}) / $DT)))
PY
)
write_mdp_md md.mdp "$NSTEPS_MD" "$DT" 300
$GMX grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 5
$GMX mdrun -deffnm md -ntmpi 1 -ntomp "$OMP_NUM_THREADS"

echo "[MD2] DONE $ROOT_DIR"
"""


def _write_executable(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def render_run_all_local(jobs: List[Job]) -> str:
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
    for j in jobs:
        sh = j.root / "run_md2.sh"
        lines.append(f'bash {sh.as_posix()}')
    return "\n".join(lines) + "\n"


def render_slurm_array(
    md_root: Path,
    jobs_list: List[Path],
    threads: int,
    slurm_partition: Optional[str],
    slurm_time: str,
    slurm_gpus: Optional[int] = None,
    slurm_mem: Optional[str] = None,
    gmxrc: Optional[str] = None,
    conda_env: Optional[str] = None,
) -> str:
    array_n = len(jobs_list)

    part_line = f"#SBATCH -p {slurm_partition}" if slurm_partition else ""
    gpu_line = f"#SBATCH --gres=gpu:{slurm_gpus}" if slurm_gpus and slurm_gpus > 0 else ""
    mem_line = f"#SBATCH --mem={slurm_mem}" if slurm_mem else ""

    # Environment bootstrap inside batch job (robust for non-interactive shells)
    
    env_lines: List[str] = []
    if gmxrc:
        env_lines.append(r'# Fix GMXRC "shell: unbound variable" under set -u')
        env_lines.append('export shell=bash')
        env_lines.append(f'if [[ -f "{gmxrc}" ]]; then source "{gmxrc}"; fi')
        env_lines.append('unset shell')
    if conda_env:
        env_lines.append(r'''
# Enable "conda activate" in non-interactive batch shell
if command -v conda >/dev/null 2>&1; then
  __conda_setup="$('conda' 'shell.bash' 'hook' 2>/dev/null)" || true
  if [[ -n "${__conda_setup:-}" ]]; then
    eval "$__conda_setup"
    conda activate ''' + conda_env + r''' || true
  else
    # Fallback: try common conda.sh location
    __conda_base="$(conda info --base 2>/dev/null || true)"
    if [[ -n "${__conda_base:-}" && -f "${__conda_base}/etc/profile.d/conda.sh" ]]; then
      source "${__conda_base}/etc/profile.d/conda.sh"
      conda activate ''' + conda_env + r''' || true
    fi
  fi
fi
'''.strip("\n"))

    env_block = "\n".join(env_lines).strip()

    # Note: jobs.list is written to md_root/jobs.list
    return f"""#!/usr/bin/env bash
#SBATCH -J mpro_md2
{part_line}
#SBATCH -t {slurm_time}
#SBATCH -c {threads}
{gpu_line}
{mem_line}
#SBATCH --array=0-{array_n-1}
#SBATCH -o {md_root.as_posix()}/slurm-%A_%a.out
#SBATCH -e {md_root.as_posix()}/slurm-%A_%a.err

set -euo pipefail
export OMP_NUM_THREADS={threads}

{env_block}

mapfile -t JOBS < { (md_root / "jobs.list").as_posix() }
ROOT="${{JOBS[$SLURM_ARRAY_TASK_ID]}}"
bash "$ROOT/run_md2.sh"
"""



def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--md_input", default="results/step6b_md_input.csv")
    ap.add_argument("--md_root", default="md")
    ap.add_argument("--prepare_only", action="store_true")
    ap.add_argument("--scheduler", choices=["local", "slurm"], default="local")
    ap.add_argument("--slurm_partition", default=None, help="If omitted, do not set #SBATCH -p (use cluster default).")
    ap.add_argument("--slurm_time", default="24:00:00")
    ap.add_argument("--threads", type=int, default=8)
    # hongmei add
    ap.add_argument("--slurm_gpus", type=int, default=None,
                help="GPUs per array task. If set, add '#SBATCH --gres=gpu:N' (or cluster policy).")
    ap.add_argument("--slurm_mem", default=None,
                    help="Memory per array task, e.g. 32G. If set, add '#SBATCH --mem=32G'.")
    ap.add_argument("--gmxrc", default=None,
                    help="Path to GMXRC to source inside the slurm job before running MD.")
    ap.add_argument("--conda_env", default=None,
                    help="Conda env name to activate inside the slurm job before running MD.")


    # QC thresholds
    ap.add_argument("--mindist_thr_nm", type=float, default=0.12)
    ap.add_argument("--fmax_thr", type=float, default=500.0)

    # fallback knobs
    ap.add_argument("--fb_emtol", type=float, default=100.0)
    ap.add_argument("--fb_em_nsteps", type=int, default=200000)
    ap.add_argument("--fb_nvt50_ps", type=float, default=20.0)
    ap.add_argument("--fb_nvt300_ps", type=float, default=80.0)

    # main equilibration lengths
    ap.add_argument("--main_nvt_ps", type=float, default=100.0)
    ap.add_argument("--npt_ps", type=float, default=100.0)

    args = ap.parse_args()

    md_input_csv = Path(args.md_input)
    md_root = Path(args.md_root)

    jobs = _read_jobs(md_input_csv, md_root)
    job_roots = [j.root for j in jobs]


    prepared: List[Job] = []
    skipped_missing_prepare = 0
    for j in jobs:
        ok = all((j.param_dir / x).exists() for x in ["complex.gro", "topol.top", "posre.itp", "ligand.itp"])
        if not ok:
            skipped_missing_prepare += 1
            continue
        prepared.append(j)

    md_root.mkdir(parents=True, exist_ok=True)

    # Write per-job runner
    for j in prepared:
        content = render_run_md2_sh(
            j,
            threads=args.threads,
            mindist_thr_nm=args.mindist_thr_nm,
            fmax_thr=args.fmax_thr,
            fb_emtol=args.fb_emtol,
            fb_em_nsteps=args.fb_em_nsteps,
            fb_nvt50_ps=args.fb_nvt50_ps,
            fb_nvt300_ps=args.fb_nvt300_ps,
            main_nvt_ps=args.main_nvt_ps,
            npt_ps=args.npt_ps,
            # conda_env=(args.conda_env.strip() if isinstance(args.conda_env, str) and args.conda_env.strip() else None),
            # gmxrc=(args.gmxrc.strip() if isinstance(args.gmxrc, str) and args.gmxrc.strip() else None),
        )
        _write_executable(j.root / "run_md2.sh", content)

    # Write master runner + jobs.list + slurm script
    _write_executable(md_root / "run_all_local_md2.sh", render_run_all_local(prepared))
    (md_root / "jobs.list").write_text("\n".join(str(j.root.resolve()) for j in prepared) + "\n", encoding="utf-8")

    slurm = render_slurm_array(
        md_root=md_root,
        jobs_list=job_roots,
        threads=args.threads,
        slurm_partition=(args.slurm_partition.strip() if isinstance(args.slurm_partition, str) and args.slurm_partition.strip() else None),
        slurm_time=args.slurm_time,
        slurm_gpus=args.slurm_gpus,
        slurm_mem=args.slurm_mem,
        gmxrc=(args.gmxrc.strip() if isinstance(args.gmxrc, str) and args.gmxrc.strip() else None),
        conda_env=(args.conda_env.strip() if isinstance(args.conda_env, str) and args.conda_env.strip() else None),
    )

    _write_executable(md_root / "slurm_array_md2.sbatch", slurm)

    print(f"[DONE] md_input={md_input_csv.resolve()}")
    print(f"[DONE] jobs prepared = {len(prepared)}  skipped_missing_prepare = {skipped_missing_prepare}")
    print(f"[DONE] local runner: { (md_root / 'run_all_local_md2.sh').resolve() }")
    print(f"[DONE] slurm script: { (md_root / 'slurm_array_md2.sbatch').resolve() }")

    if args.prepare_only:
        return 0

    if args.scheduler == "local":
        os.execv("/bin/bash", ["bash", str((md_root / "run_all_local_md2.sh").resolve())])
    else:
        print("[INFO] Submit with: sbatch md/slurm_array_md2.sbatch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
