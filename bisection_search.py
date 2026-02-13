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

# Bisection search parameters
BOUND_TYPE = 'upper'  # 'upper' or 'lower'
INITIAL_LOWER = 0.0
INITIAL_UPPER = 1.0
TOLERANCE = 0.01
MAX_ITERATIONS = 20

# Path to executable
EXECUTABLE_PATH = './build/implicit_learning'

#############################################


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
        
        print(f"  Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            if result.returncode != 0:
                print(f"  ERROR: C++ program failed with exit code {result.returncode}")
                print(f"  stderr: {result.stderr}")
                return False
            return True
        except subprocess.TimeoutExpired:
            print(f"  ERROR: C++ program timed out after 1800 seconds")
            return False
        except Exception as e:
            print(f"  ERROR: Failed to run C++ program: {e}")
            return False
    
    def check_feasibility(self):
        # STUB: Will be replaced with actual solver call
        print("  [STUB] Checking feasibility... (always returning True for now)")
        return True
    
    def search(self):

        print(f"\n{'='*60}")
        print(f"Starting Bisection Search")
        print(f"{'='*60}")
        print(f"Relation: {self.relation_name}")
        print(f"Bound Type: {self.bound_type}")
        print(f"Initial Range: [{self.lower}, {self.upper}]")
        print(f"Tolerance: {self.tolerance}")
        print(f"Fixed Gene: {self.fixed_gene}")
        print(f"Fixed Enzyme: {self.fixed_enzyme}")
        print(f"Data File: {self.data_file}")
        print(f"{'='*60}\n")
        
        iteration = 0
        
        while (self.upper - self.lower) > self.tolerance and iteration < self.max_iterations:
            iteration += 1
            mid = (self.lower + self.upper) / 2.0
            
            print(f"Iteration {iteration}/{self.max_iterations}:")
            print(f"  Current range: [{self.lower:.6f}, {self.upper:.6f}]")
            print(f"  Testing midpoint: {mid:.6f}")
            
            # Run C++ program with current bound
            if not self.run_cpp_program(mid):
                print(f"  C++ program failed, aborting search")
                return None
            
            # Check if solution is feasible
            is_feasible = self.check_feasibility()
            
            # Update bounds based on feasibility
            if self.bound_type == 'upper':
                if is_feasible:
                    # Can tighten upper bound
                    self.upper = mid
                    print(f"  ✓ Feasible - tightening upper bound to {mid:.6f}")
                else:
                    # Need to relax upper bound
                    self.lower = mid
                    print(f"  ✗ Infeasible - relaxing upper bound, new lower={mid:.6f}")
            else:  # lower bound
                if is_feasible:
                    # Can tighten lower bound
                    self.lower = mid
                    print(f"  ✓ Feasible - tightening lower bound to {mid:.6f}")
                else:
                    # Need to relax lower bound
                    self.upper = mid
                    print(f"  ✗ Infeasible - relaxing lower bound, new upper={mid:.6f}")
            
            print()
        
        # Return final bound
        final_bound = self.upper if self.bound_type == 'upper' else self.lower
        
        print(f"{'='*60}")
        print(f"Search Complete!")
        print(f"{'='*60}")
        print(f"Final {self.bound_type} bound: {final_bound:.6f}")
        print(f"Iterations: {iteration}")
        print(f"Final range: [{self.lower:.6f}, {self.upper:.6f}]")
        print(f"{'='*60}\n")
        
        return final_bound


def main():
    
    print("="*60)
    print("CONFIGURATION (edit at top of script)")
    print("="*60)
    print(f"Fixed Gene:     {FIXED_GENE}")
    print(f"Fixed Enzyme:   {FIXED_ENZYME}")
    print(f"Data File:      {DATA_FILE}")
    print(f"Relation:       {RELATION_NAME}")
    print(f"Bound Type:     {BOUND_TYPE}")
    print(f"Initial Range:  [{INITIAL_LOWER}, {INITIAL_UPPER}]")
    print(f"Tolerance:      {TOLERANCE}")
    print(f"Max Iterations: {MAX_ITERATIONS}")
    print("="*60)
    
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