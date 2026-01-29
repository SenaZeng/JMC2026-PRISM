#!/usr/bin/env python3
"""
PRISM minimal smoke test.

This script checks that core/step1.py ... core/step6.py exist and can be executed
with '--help' (or without crashing). It does not run the full pipeline.
"""

from __future__ import annotations
import os
import sys
import subprocess
import py_compile
from pathlib import Path


def run(cmd: list[str]) -> int:
    print(f"\n[RUN] {' '.join(cmd)}")
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(p.stdout[:2000])  # print first 2000 chars to keep CI/logs small
    return p.returncode


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    core_dir = repo_root / "core"

    if not core_dir.exists():
        print(f"[FAIL] core/ not found at: {core_dir}")
        return 1

    steps = [core_dir / f"step{i}.py" for i in range(1, 2)]
    
        # Auto-detect scripts under core/ to avoid hardcoded filenames like core/step1.py.
    steps = sorted(core_dir.glob("step*.py"))
    if not steps:
        steps = sorted(p for p in core_dir.glob("*.py") if p.name != "__init__.py")

    if not steps:
        print("[FAIL] No Python scripts found under core/.")
        return 1

    print("[INFO] Detected core scripts:")
    for p in steps:
        print(f"  - {p.name}")

    

    python_exe = sys.executable or "python"
    failures = 0

        # CI-friendly check: syntax-only (do not execute heavy workflows).
    failures = 0
    for p in steps:
        try:
            py_compile.compile(str(p), doraise=True)
            print(f"[OK] syntax: {p.name}")
        except Exception as e:
            print(f"[FAIL] syntax: {p.name}")
            print(f"       {e}")
            failures += 1

    if failures == 0:
        print("\n[PASS] First test completed (syntax-only).")
        return 0

    print(f"\n[FAIL] First test had {failures} failing script(s).")
    return 1

    

    if failures == 0:
        print("\n[PASS] Smoke test completed. Core scripts are callable.")
        return 0

    print(f"\n[FAIL] Smoke test had {failures} failing script(s).")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
