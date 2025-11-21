#include "executor.h"
#include "config.h"
#include "domain.h"
#include "metrics.h"

#include <iostream>
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