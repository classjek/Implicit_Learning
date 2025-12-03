#include "executor.h"
#include "config.h"
#include "domain.h"
#include "metrics.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo> // for debugging, can remove later

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
  // Convert to vector of strings
  std::vector<std::vector<std::string>> typedGroundNames;
  typedGroundNames.push_back(std::vector<std::string>(groundNames.genes.begin(), groundNames.genes.end()));
  typedGroundNames.push_back(std::vector<std::string>(groundNames.enzymes.begin(), groundNames.enzymes.end()));
  typedGroundNames.push_back(std::vector<std::string>(groundNames.reactions.begin(), groundNames.reactions.end()));
  typedGroundNames.push_back(std::vector<std::string>(groundNames.compounds.begin(), groundNames.compounds.end()));

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
  std::cout << '\n' << "Printing All Constraints(" << universal_constraints.size() << "):" << std::endl;
  for (size_t i = 0; i < universal_constraints.size(); i++) {
      std::cout << "    constraint[" << i << "] takes " << "4" << " args: "  << universal_constraints[i].poly.toString() << " " << (universal_constraints[i].cmp == kb::Cmp::GE0 ? ">=" : "=") << " 0\n";
  }
  

std::unordered_map<size_t,int> groundMap; 
std::vector<std::vector<std::vector<int>>> finalResults(universal_constraints.size());

// Build smaller set of groundNames for testing
std::vector<std::vector<std::string>> groundNamesTest(typedGroundNames.size());
groundNamesTest[0].assign(typedGroundNames[0].begin(), typedGroundNames[0].begin()+400);
groundNamesTest[1].assign(typedGroundNames[1].begin(), typedGroundNames[1].begin()+400);
groundNamesTest[2].assign(typedGroundNames[2].begin(), typedGroundNames[2].begin()+400);
groundNamesTest[3].assign(typedGroundNames[3].begin(), typedGroundNames[3].begin()+200);
for (auto& elem : groundNamesTest) { std::cout << elem.size() << ", "; }
std::cout << std::endl;

// to avoid storing partially ground results, use nested DFS
// consider one constraint at a time, generate all partial observations for it
// So we have a function that, given a constraint, returns what it takes as input 
// domain::generateGrounding(universal_constraints, typedGroundNames, groundMap, finalResults);
domain::generateGrounding(universal_constraints, groundNamesTest, groundMap, finalResults); // for Testing

cp.tick("After grounding"); 

std::cout << "Grounded Atom Map (total " << groundMap.size() << " atoms):" << std::endl;
std::cout << "finalResults size: " << finalResults.size() << std::endl;
std::cout << "finalResults[0] size: " << finalResults[0].size() << std::endl;
std::cout << "finalResults[1] size: " << finalResults[1].size() << std::endl;

// discard groundMap
groundMap.clear();
groundMap = std::unordered_map<size_t, int>();


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
