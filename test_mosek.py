#!/usr/bin/env python3

import subprocess
import sys
import os
import csv
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, TimeoutError as FuturesTimeout

from run_mosek import solve as mosek_solve

DATA_FILE       = '1483249_new.pl'
EXECUTABLE_PATH = './implicit_learning'
DAT_FILE        = Path('./data/sparsepop_output_test.dat-s')
TSV_FILE        = './data/jobArray_function.tsv'

BOUND_TYPE      = 'upper'
BOUND_VALUE     = 0.01

MOSEK_TIMEOUT   = 600   
CPP_TIMEOUT     = 300  

OUTPUT_CSV      = 'probe_results.csv'


def load_all_jobs(tsv_file):
    jobs = []
    with open(tsv_file) as f:
        lines = f.readlines()
    for i, line in enumerate(lines[1:], start=1):   # skip header
        parts = line.strip().split('\t')
        if len(parts) != 3:
            continue
        fixed_gene    = 'g' + parts[0]
        fixed_enzyme  = 'ec_' + parts[1].replace('.', '_').replace('-', '')
        relation_name = parts[2]
        jobs.append((i, fixed_gene, fixed_enzyme, relation_name))
    return jobs


def run_cpp(relation_name, fixed_gene, fixed_enzyme):
    cmd = [
        EXECUTABLE_PATH,
        '--bound-atom',  relation_name,
        '--bound-value', str(BOUND_VALUE),
        '--bound-type',  BOUND_TYPE,
        '--fixedGene',   fixed_gene,
        '--fixedEnzyme', fixed_enzyme,
        '--fileName',    DATA_FILE,
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=CPP_TIMEOUT, cwd='build'
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        return False
    except Exception:
        return False


def _mosek_worker(dat_path):
    """Runs in a separate process so we can kill it on timeout."""
    prob = mosek_solve(dat_path, verbose=False)
    return prob.status


def check_mosek_with_timeout(dat_file):
    """Returns (status_str, elapsed_seconds).
    status_str: 'optimal' | 'infeasible' | 'other:<status>' | 'TIMEOUT' | 'ERROR'
    """
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=1) as ex:
        future = ex.submit(_mosek_worker, str(dat_file.absolute()))
        try:
            status = future.result(timeout=MOSEK_TIMEOUT)
            elapsed = time.perf_counter() - t0
            if status == 'optimal':
                return 'OK_feasible', elapsed
            elif 'infeasible' in status:
                return 'OK_infeasible', elapsed
            else:
                return f'OK_other:{status}', elapsed
        except FuturesTimeout:
            elapsed = time.perf_counter() - t0
            return 'TIMEOUT', elapsed
        except Exception as e:
            elapsed = time.perf_counter() - t0
            return f'ERROR:{e}', elapsed


def main():
    if not os.path.exists(TSV_FILE):
        print(f"TSV file not found: {TSV_FILE}")
        sys.exit(1)

    jobs = load_all_jobs(TSV_FILE)
    print(f"Loaded {len(jobs)} jobs from {TSV_FILE}")
    print(f"Bound: {BOUND_TYPE} {BOUND_VALUE}  |  MOSEK timeout: {MOSEK_TIMEOUT}s\n")

    with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['job', 'relation', 'fixed_gene', 'fixed_enzyme', 'cpp_ok', 'mosek_status', 'mosek_secs'])

        for (job_num, fixed_gene, fixed_enzyme, relation_name) in jobs:
            print(f"[{job_num}/{len(jobs)}] {relation_name} ...", end=' ', flush=True)

            # Run grounding
            cpp_ok = run_cpp(relation_name, fixed_gene, fixed_enzyme)
            if not cpp_ok:
                print("CPP_FAIL")
                writer.writerow([job_num, relation_name, fixed_gene, fixed_enzyme, False, 'CPP_FAIL', ''])
                csvfile.flush()
                continue

            # Check .dat-s file existence
            if not DAT_FILE.exists():
                print("NO_FILE")
                writer.writerow([job_num, relation_name, fixed_gene, fixed_enzyme, True, 'NO_FILE', ''])
                csvfile.flush()
                continue

            # Attempt MOSEK solve
            mosek_status, elapsed = check_mosek_with_timeout(DAT_FILE)
            print(f"{mosek_status}  ({elapsed:.1f}s)")
            writer.writerow([job_num, relation_name, fixed_gene, fixed_enzyme, True, mosek_status, f'{elapsed:.1f}'])
            csvfile.flush()

    print(f"\nDone. Results written to {OUTPUT_CSV}")

if __name__ == '__main__':
    main()