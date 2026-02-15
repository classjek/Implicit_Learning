#!/usr/bin/env python3

import subprocess
import sys
import os
from pathlib import Path

## HARDCODED CONFIG FOR TESTING

# Fixed parameters for grounding
FIXED_GENE = "g100036608"
FIXED_ENZYME = "ec_3_1_3_48"
DATA_FILE = "R-HSA-1483249_data.pl"

# Relation to find bounds for
RELATION_NAME = "function(g100036608,ec_3_4_21)"

LOG_FILE = "bisection_search" + RELATION_NAME + ".log";

# Bisection search parameters
BOUND_TYPE = 'upper'  # 'upper' or 'lower'
INITIAL_LOWER = 0.0
INITIAL_UPPER = 1.0
TOLERANCE = 0.01
MAX_ITERATIONS = 100

# Path to executable
EXECUTABLE_PATH = './implicit_learning'

#############################################

def log_print(message="", end="\n"):
    """Print to both console and log file."""
    print(message, end=end)
    with open(LOG_FILE, 'a') as f:
        f.write(message + end)


class BisectionSearch:
    def __init__(self, relation_name, bound_type, initial_lower, initial_upper, tolerance, max_iterations,
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
    
    log_print("="*60)
    log_print("CONFIGURATION")
    log_print("="*60)
    log_print(f"Fixed Gene:     {FIXED_GENE}")
    log_print(f"Fixed Enzyme:   {FIXED_ENZYME}")
    log_print(f"Data File:      {DATA_FILE}")
    log_print(f"Relation:       {RELATION_NAME}")
    log_print(f"Bound Type:     {BOUND_TYPE}")
    log_print(f"Initial Range:  [{INITIAL_LOWER}, {INITIAL_UPPER}]")
    log_print(f"Tolerance:      {TOLERANCE}")
    log_print(f"Max Iterations: {MAX_ITERATIONS}")
    log_print("="*60)
    
    # Create and run search
    search = BisectionSearch(
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
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()