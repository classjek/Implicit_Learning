#!/usr/bin/env python3

import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

## HARDCODED CONFIG FOR TESTING
# Fixed parameters for grounding
# FIXED_GENE = "g100036608"
# FIXED_ENZYME = "ec_3_1_3_48"
# RELATION_NAME = "function(g100036608,ec_3_4_21)"
# LOG_FILE = "bisection_search" + RELATION_NAME + ".log"

DATA_FILE = "R-HSA-1483249_data.pl"

# Bisection search parameters
BOUND_TYPE = 'upper'  # 'upper' or 'lower'
INITIAL_LOWER = 0.0
INITIAL_UPPER = 1.0
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# Path to executable
EXECUTABLE_PATH = './implicit_learning'
#############################################

def load_job_config(job_number, tsv_file='./data/jobArray_function.tsv'):
    if not os.path.exists(tsv_file):
        raise FileNotFoundError(f"Job array file not found: {tsv_file}")
    
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
    
    data_index = job_number  # job_number is 1-indexed, matches line 1 after header
    
    if data_index < 1 or data_index >= len(lines):
        raise ValueError(f"Invalid job number {job_number}. Valid range: 1-{len(lines)-1}")
    
    # Parse TSV line
    parts = lines[data_index].strip().split('\t')
    if len(parts) != 3:
        raise ValueError(f"Invalid TSV format at line {data_index}: expected 3 columns")
    
    # Format the fields
    fixed_gene = 'g' + parts[0]  # Add 'g' prefix
    fixed_enzyme = 'ec_' + parts[1].replace('.', '_').replace('-', '') 
    ground_relation = parts[2]
    
    return fixed_gene, fixed_enzyme, ground_relation

LOG_FILE = None  # Will be set in main() based on job number

# print to both console and log file
def log_print(message="", end="\n"):
    print(message, end=end)
    if LOG_FILE is not None:
        with open(LOG_FILE, 'a') as f:
            f.write(message + end)


class BisectionSearch:
    def __init__(self, job_number, relation_name, bound_type, initial_lower, initial_upper, tolerance, max_iterations,
                 fixed_gene, fixed_enzyme, data_file, executable_path):
        """
            relation_name: Name of the relation (e.g., "function(g100036608,ec_3_4_21)")
            bound_type: 'upper' or 'lower'
            initial_lower: Starting lower bound
            initial_upper: Starting upper bound
            tolerance: Convergence threshold
            max_iterations: Maximum number of iterations
            fixed_gene: Fixed gene for grounding
            fixed_enzyme: Fixed enzyme for grounding
            data_file: Data file name
            executable_path: Path to compiled C++ executable
        """
        self.job_number = job_number
        self.relation_name = relation_name
        self.bound_type = bound_type
        self.lower = initial_lower
        self.upper = initial_upper
        self.tolerance = tolerance
        self.max_iterations = max_iterations
        self.fixed_gene = fixed_gene
        self.fixed_enzyme = fixed_enzyme
        self.data_file = data_file
        self.executable_path = executable_path
        self.log_file = f"bisection_search_job{job_number}.log"  # Job-specific log
        
    # Run with a specified bound
    def run_cpp_program(self, bound_value):
        cmd = [
            self.executable_path,
            '--bound-atom', self.relation_name,
            '--bound-value', str(bound_value),
            '--bound-type', self.bound_type,
            '--fixedGene', self.fixed_gene,
            '--fixedEnzyme', self.fixed_enzyme,
            '--fileName', self.data_file
        ]
        
        log_print(f"  Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800, cwd='build')
            if result.returncode != 0:
                log_print(f"  ERROR: C++ program failed with exit code {result.returncode}")
                log_print(f"  stderr: {result.stderr}")
                return False
            return True
        except subprocess.TimeoutExpired:
            log_print(f"  ERROR: C++ program timed out after 1800 seconds")
            return False
        except Exception as e:
            log_print(f"  ERROR: Failed to run C++ program: {e}")
            return False
    
    # Todo
    def check_feasibility(self, culorads_timeout=3600):
        output_file = Path('./data/sparsepop_output_test.dat-s')

        if not output_file.exists():
            log_print(f"  Output file not found: {output_file}")
            return False
    
        file_size = output_file.stat().st_size
        log_print(f"  Found output file: {output_file} ({file_size} bytes)")
    
        # Run cuLoRADS
        culorads_cmd = [
            './bin/cuLoRADS',
            '--filePath', str(output_file.absolute()),
            '--timeSecLimit', str(culorads_timeout)
        ]
        log_print(f"  Running cuLoRADS (timeout: {culorads_timeout}s)...")
    
        try:
            result = subprocess.run(
                culorads_cmd,
                capture_output=True,
                text=True,
                timeout=culorads_timeout + 10
            )

            # Check if cuLoRADS succeeded
            if result.returncode != 0:
                log_print(f"cuLoRADS failed with exit code {result.returncode}")
                log_print(f"  stderr: {result.stderr[:200]}")  # First 200 chars
                output_file.unlink()
                return False

            # Parse output to determine feasibility
            output = result.stdout

            # Look for the success marker
            if "Problem Solved" in output:
                log_print(f"cuLoRADS solved successfully - FEASIBLE")
                is_feasible = True
            else:
                log_print(f"cuLoRADS did not solve - INFEASIBLE")
                is_feasible = False

            # Cleanup
            output_file.unlink()
            log_print(f"Deleted output file")
        
            return is_feasible
        
        except subprocess.TimeoutExpired:
            log_print(f"cuLoRADS timed out after {culorads_timeout}s - treating as INFEASIBLE")
            output_file.unlink()
            return False
        except Exception as e:
            log_print(f"  Error running cuLoRADS: {e}")
            output_file.unlink()
            return False
    
    def search(self):

        log_print(f"\n{'='*60}")
        log_print(f"Starting Bisection Search")
        log_print(f"{'='*60}")
        log_print(f"Relation: {self.relation_name}")
        log_print(f"Bound Type: {self.bound_type}")
        log_print(f"Initial Range: [{self.lower}, {self.upper}]")
        log_print(f"Tolerance: {self.tolerance}")
        log_print(f"Fixed Gene: {self.fixed_gene}")
        log_print(f"Fixed Enzyme: {self.fixed_enzyme}")
        log_print(f"Data File: {self.data_file}")
        log_print(f"{'='*60}\n")
        
        iteration = 0
        
        while (self.upper - self.lower) > self.tolerance and iteration < self.max_iterations:
            iteration += 1
            mid = (self.lower + self.upper) / 2.0
            
            log_print(f"Iteration {iteration}/{self.max_iterations}:")
            log_print(f"  Current range: [{self.lower:.6f}, {self.upper:.6f}]")
            log_print(f"  Testing midpoint: {mid:.6f}")
            
            # Run C++ program with current bound
            if not self.run_cpp_program(mid):
                log_print(f"  C++ program failed, aborting search")
                return None
            
            # Check if solution is feasible
            is_feasible = self.check_feasibility()
            
            # Update bounds based on feasibility
            if self.bound_type == 'upper':
                if is_feasible:
                    # Can tighten upper bound
                    self.upper = mid
                    log_print(f"  ✓ Feasible - tightening upper bound to {mid:.6f}")
                else:
                    # Need to relax upper bound
                    self.lower = mid
                    log_print(f"  ✗ Infeasible - relaxing upper bound, new lower={mid:.6f}")
            else:  # lower bound
                if is_feasible:
                    # Can tighten lower bound
                    self.lower = mid
                    log_print(f"  ✓ Feasible - tightening lower bound to {mid:.6f}")
                else:
                    # Need to relax lower bound
                    self.upper = mid
                    log_print(f"  ✗ Infeasible - relaxing lower bound, new upper={mid:.6f}")
            
            log_print()
        
        # Return final bound
        final_bound = self.upper if self.bound_type == 'upper' else self.lower
        
        log_print(f"{'='*60}")
        log_print(f"Search Complete!")
        log_print(f"{'='*60}")
        log_print(f"Final {self.bound_type} bound: {final_bound:.6f}")
        log_print(f"Iterations: {iteration}")
        log_print(f"Final range: [{self.lower:.6f}, {self.upper:.6f}]")
        log_print(f"{'='*60}\n")
        
        return final_bound


def main():    
    global LOG_FILE 

    # Parse command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python3 bisection_search.py <job_number>")
        print("\nExample:")
        print("  python3 bisection_search.py 1")
        print("\nJob configurations are read from jobArray_function.tsv")
        sys.exit(1)
    
    try:
        job_number = int(sys.argv[1])
    except ValueError:
        print(f"ERROR: Job number must be an integer, got: {sys.argv[1]}")
        sys.exit(1)
    
    # Load job-specific configuration from TSV
    try:
        FIXED_GENE, FIXED_ENZYME, RELATION_NAME = load_job_config(job_number)
    except Exception as e:
        print(f"ERROR loading job config: {e}")
        sys.exit(1)

    LOG_FILE = f"bisection_search_job{job_number}.log"

    with open(LOG_FILE, 'w') as f:
        f.write("="*60 + "\n")
        f.write(f"BISECTION SEARCH LOG - Job {job_number}\n")
        f.write(f"Started: {datetime.now()}\n")
        f.write("="*60 + "\n\n")
    
    print("="*60)
    print(f"JOB {job_number} CONFIGURATION")
    print("="*60)
    print(f"Fixed Gene:     {FIXED_GENE}")
    print(f"Fixed Enzyme:   {FIXED_ENZYME}")
    print(f"Data File:      {DATA_FILE}")
    print(f"Relation:       {RELATION_NAME}")
    print(f"Bound Type:     {BOUND_TYPE}")
    print(f"Initial Range:  [{INITIAL_LOWER}, {INITIAL_UPPER}]")
    print(f"Tolerance:      {TOLERANCE}")
    print(f"Max Iterations: {MAX_ITERATIONS}")
    print(f"Log File:       {LOG_FILE}")
    print("="*60)
    
    # Create and run search
    search = BisectionSearch(
        job_number=job_number,
        relation_name=RELATION_NAME,
        bound_type=BOUND_TYPE,
        initial_lower=INITIAL_LOWER,
        initial_upper=INITIAL_UPPER,
        tolerance=TOLERANCE,
        max_iterations=MAX_ITERATIONS,
        fixed_gene=FIXED_GENE,
        fixed_enzyme=FIXED_ENZYME,
        data_file=DATA_FILE,
        executable_path=EXECUTABLE_PATH
    )
    
    final_bound = search.search()
    
    if final_bound is not None:
        # Write final result to a summary file
        with open(f"result_job{job_number}.txt", 'w') as f:
            f.write(f"Job: {job_number}\n")
            f.write(f"Relation: {RELATION_NAME}\n")
            f.write(f"Fixed Gene: {FIXED_GENE}\n")
            f.write(f"Fixed Enzyme: {FIXED_ENZYME}\n")
            f.write(f"Final {BOUND_TYPE} bound: {final_bound:.6f}\n")
        print(f"\n✓ Result written to result_job{job_number}.txt")
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()