#include "executor.h"
#include "config.h"
#include "domain.h"
#include "metrics.h"

#include <iostream>
#include <fstream>
#include <cmath>

// int main(int argc, char** argv) {
int main() {
  metrics::Checkpoint cp("Program Start"); // For metrics tracking

  Config cfg;  // keep defaults for now

  // Load in Observations 
  // store them in some list or vector
  // consider all constraints together, add generics and generate equivalence classes
  // consider generating equivalence classes in parallel
  // Generate maps from each constraint -> equivalence class
  domain::initializePredicateSignatures(); 
  domain::GroundNames groundNames;

  // createProbLog Parser
  domain::ProbLogParser parser(groundNames); 
  // Open file and parse constraints
  std::string filename = "../data/groundFacts.pl";
  std::vector<kb::Constraint> constraints = parser.parseFile(filename);
  cp.tick("After parsing"); 

  // at this point, we should have our groundNames structs populated and our constraints vector filled 
  std::cout << "Parsed " << constraints.size() << " constraints" << std::endl;
  std::cout << "Ground Genes: " << groundNames.genes.size() << " , Enzymes: " << groundNames.enzymes.size() 
            << " , Reactions: " << groundNames.reactions.size() << " , Compounds: " << groundNames.compounds.size() << std::endl;

  // Read in Universally Quantified Constraints
  std::ifstream in("../data/universalConstraints.txt");   
  if (!in) {
      std::cerr << "Could not open constraints file\n";
      return 1;
  }

  std::vector<kb::Constraint> universal_constraints;
  std::string line;
  while (std::getline(in, line)) {
        if (line.empty() || (line[0] == '#' && line[1] == '#')) continue;    
        try {
            // Only add constraint if not already present
            auto constraint = domain::parseConstraint(line);
            if (std::find(universal_constraints.begin(), universal_constraints.end(), constraint) == universal_constraints.end()) {
                universal_constraints.push_back(constraint);
            }
        } catch (const std::exception& e) {
            std::cerr << "Parse error in line: \"" << line << "\"\n  " << e.what() << '\n';
        }
  }
  in.close();
  // std::cout << "Universal constraints added: " << universal_constraints.size() << std::endl;
  std::cout << '\n' << "Printing All Constraints(" << constraints.size() << "):" << std::endl;
  for (size_t i = 0; i < universal_constraints.size(); i++) {
      std::cout << "    constraint[" << i << "] takes " << "6" << " args: "  << universal_constraints[i].poly.toString() << " " << (universal_constraints[i].cmp == kb::Cmp::GE0 ? ">=" : "=") << " 0\n";
  }


  int num_observations = 6;
  cfg.omp_threads = std::floor(36 / num_observations);

  std::cout << "\nThreads: " << cfg.omp_threads << std::endl;
  // Later pass vector of partialobservations and equivalence class map to executor
  Executor ex(cfg);

  // run should return information about bounds of each equivalence class per partial observation
  // Depends how equivalence classes are stored
  int rc = ex.run();
  std::cout << "[main] done, rc=" << rc << "\n" << std::endl;

  cp.tick("Program End"); 
  cp.print(); 

  return rc;
}