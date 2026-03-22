#!/usr/bin/env python3
"""
test_bound.py — manually test a single bound constraint.

Usage:
  python3 test_bound.py 'function(g504545,g100036608)' upper 0.5
  python3 test_bound.py 'function(g504545,g100036608)' lower 0.3
"""

import subprocess
import sys
from pathlib import Path
from run_mosek import solve as mosek_solve

# ── Config (edit these to match your setup) ──────────────────────
DATA_FILE       = '1483249_new.pl'
# FIXED_GENE      = 'g100036608'
# FIXED_GENE = 'g100038240'
# FIXED_ENZYME    = 'ec_1_13_99_1'
FIXED_GENE      = 'g10295458242'
FIXED_ENZYME = 'ec_2_7_1_134'
# FIXED_ENZYME    = 'ec_3_6_1_61'
EXECUTABLE_PATH = './implicit_learning'
DAT_FILE        = Path('./data/sparsepop_output_test.dat-s')
# ─────────────────────────────────────────────────────────────────


def run_cpp(atom_name, bound_type, bound_value):
    cmd = [
        EXECUTABLE_PATH,
        '--bound-atom',   atom_name,
        '--bound-value',  str(bound_value),
        '--bound-type',   bound_type,
        '--fixedGene',    FIXED_GENE,
        '--fixedEnzyme',  FIXED_ENZYME,
        '--fileName',     DATA_FILE,
    ]
    print(f"\n[CPP] Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False, text=True, timeout=1800, cwd='build')
    if result.returncode != 0:
        print(f"[CPP] ERROR: exit code {result.returncode}")
        return False
    return True


def check_feasibility():
    if not DAT_FILE.exists():
        print(f"[MOSEK] dat-s file not found: {DAT_FILE}")
        return None
    print(f"\n[MOSEK] Solving {DAT_FILE} ({DAT_FILE.stat().st_size} bytes) ...")
    prob = mosek_solve(str(DAT_FILE.absolute()), verbose=True)
    print(f"\n[RESULT] status={prob.status}  obj={prob.value}")
    return prob.status == 'optimal'

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 test_bound.py <atom_name> <upper|lower> <value>")
        print("Example: python3 test_bound.py 'function(g504545,g100036608)' upper 0.5")
        sys.exit(1)

    atom_name   = sys.argv[1]
    bound_type  = sys.argv[2]
    bound_value = float(sys.argv[3])

    print(f"\n{'='*55}")
    print(f"Testing:  {atom_name} {bound_type} {bound_value}")
    print(f"{'='*55}")

    ok = run_cpp(atom_name, bound_type, bound_value)
    if not ok:
        print("[RESULT] C++ failed — cannot determine feasibility")
        sys.exit(1)

    feasible = check_feasibility()
    print(f"\n{'='*55}")
    if feasible is True:
        print(f"  ✓  FEASIBLE   — {atom_name} {bound_type} bound {bound_value} is consistent")
    elif feasible is False:
        print(f"  ✗  INFEASIBLE — {atom_name} {bound_type} bound {bound_value} is violated")
    else:
        print(f"  ?  UNKNOWN    — solver did not return a clean status")
    print(f"{'='*55}\n")