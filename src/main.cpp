#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo> // for debugging, can remove later

#include "config.h"
#include "domain.h"
#include "executor.h"
#include "metrics.h"
#include "spop.h"
#include "streaming.h"

// int main(int argc, char** argv) {
int main(int argc, char** argv) {
    metrics::Checkpoint cp("Program Start"); // For metrics tracking

    // Example cmd line: 
    //        hom./implicit_learning --bound-atom "function(g100036608,ec_3_4_21)" --bound-value 0.75 --bound-type upper --fixedGene "g100036608" --fixedEnzyme "ec_3_1_3_48" --fileName "R-HSA-1483249_data.pl"

    // Parse command-line arguments for bound constraint
    std::string DATA_FILE = ""; 
    std::string fixedGene = ""; std::string fixedEnzyme = ""; 
    std::string cl_atomName = ""; double boundValue = 0.5;
    bool isLower = true;  // true for >=, false for <=
    std::vector<domain::BoundConstraint> cl_bounds;

    // Simple command-line parsing
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--bound-atom" && i + 1 < argc) {
            cl_atomName = argv[++i];
        } else if (arg == "--bound-value" && i + 1 < argc) {
            boundValue = std::atof(argv[++i]);
        } else if (arg == "--bound-type" && i + 1 < argc) {
            std::string type = argv[++i];
            isLower = (type == "lower"); 
        } else if (arg == "--fixedGene" && i + 1 < argc) {
            fixedGene = argv[++i];
        } else if (arg == "--fixedEnzyme" && i + 1 < argc) {
            fixedEnzyme = argv[++i];
        } else if (arg == "--fileName" && i + 1 < argc) {
            DATA_FILE = argv[++i];
        }
    }   

    std::cout << "Received cmd line arguments:" << std::endl;
    std::cout << "  Fixed Gene: " << fixedGene << std::endl;
    std::cout << "  Fixed Enzy: " << fixedEnzyme << std::endl;
    std::cout << "  Atom name: " << cl_atomName << std::endl;
    std::cout << "  Bound value: " << boundValue << std::endl;
    std::cout << "  Bound type: " << (isLower ? "lower (>=)" : "upper (<=)") << std::endl;


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
  //   std::string filename = "../data/groundFacts.pl";
  if (DATA_FILE == "") { 
    std::cout << "Warning: No data file indicated. Using default." << std::endl;
    DATA_FILE = "1483249_new.pl";
  }
  std::string filename = "../data/" + DATA_FILE;

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

  if (fixedGene == "" && fixedEnzyme == ""){
    std::cout << " - Warning: No fixedGene and fixedEnzyme specified, using default" << std::endl;
    fixedGene = "g100036608"; fixedEnzyme = "g100036608"; 
  }
  typedGroundNames.push_back(std::vector<std::string>{fixedGene}); 
  typedGroundNames.push_back(std::vector<std::string>{fixedEnzyme}); 
  if (cl_atomName == "") { 
    std::cout << " - Warning: No bounded atom specified, using default" << std::endl;
    cl_atomName == "function(g100036608,ec_3_4_21)";
  }

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
groundNamesTest[0].assign(typedGroundNames[0].begin(), typedGroundNames[0].begin()+200); // genes 100
groundNamesTest[1].assign(typedGroundNames[1].begin(), typedGroundNames[1].begin()+25); // enzymes 27
groundNamesTest[2].assign(typedGroundNames[2].begin(), typedGroundNames[2].begin()+1); // reactions
groundNamesTest[3].assign(typedGroundNames[3].begin(), typedGroundNames[3].begin()+1); //compounds

groundNamesTest[0].push_back("g100036608");  
groundNamesTest[0].push_back("g100037840");  
groundNamesTest[1].push_back("ec_3_1_3_48"); 
groundNamesTest[1].push_back("ec_2_3_2");     
// groundNamesTest[3].push_back("ec_2_7_1_134");


std::cout << std::endl;

// for (auto& elem : groundNamesTest[0]){
//     std::cout << " - [gene] " << elem << std::endl;
// }
// for (auto& elem : groundNamesTest[1]){
//     std::cout << " - [enzyme] " << elem << std::endl;
// }

for (auto& elem : groundNamesTest) { std::cout << elem.size() << ", "; }

std::cout << std::endl;

// ground universally quantified constraints
// domain::generateGrounding(universal_constraints, typedGroundNames, groundMap, finalResults); // for Testing
domain::generateGrounding(universal_constraints, groundNamesTest, groundMap, finalResults); // for Testing
cp.tick("After grounding"); 

std::vector<domain::BoundConstraint> bounds; 
size_t hash = std::hash<std::string>{}(cl_atomName);
auto it = groundMap.find(hash);
if (it == groundMap.end()){
    std::cerr << "Warning! Trying to place bound on unknown ground atom." << std::endl;
} else {
    bounds.push_back({
        it->second, // atomID
        boundValue,
        isLower 
    });
    std::string boundType = isLower ? ">=" : "<=";
    std::cout << " - Added bound: " << cl_atomName << " " << boundType << " " << boundValue << " (atomID=" << it->second << ") - \n" << std::endl;
}

// build observed values from facts
// Build observed values from ground facts
std::cout << "We have " << constraints.size() << " constraints" << std::endl;
std::vector<double> observedValueById = domain::buildObservedValues(constraints, groundMap, groundMap.size());

std::vector<int> polyWidth; // holds the number of arguments taken by polynomial i
std::vector<int> gndOff; // holds offset used to access the gndData for each polynomial
std::vector<int> gndData; // every valid grounding vector, stored contiguously
domain::createGroundingRepresentation(finalResults, polyWidth, gndOff, gndData);
cp.tick("After Sparse Rep"); 

std::cout << "Grounded Atom Map (total " << groundMap.size() << " atoms):" << std::endl;
std::cout << "finalResults size: " << finalResults.size() << std::endl;
std::cout << "finalResults[0] size: " << finalResults[0].size() << std::endl;
std::cout << "finalResults[1] size: " << finalResults[1].size() << "\n" <<std::endl;

int newNumVars = groundMap.size(); // number of new variables 
int newNumConst = 0; // number of grounded constraints
for (const auto& constr : finalResults) { newNumConst += constr.size(); }

// discard groundMap
groundMap.clear();
groundMap = std::unordered_map<size_t, int>();
// discard finalResults
finalResults.clear(); 
finalResults.shrink_to_fit();
cp.tick("After clearing");

// extract constraint format from grounded constraints and include that here
// I believe writeGMS fills in with bounds that will all be replaced later, TODO: confirm
std::string fileName = domain::writeGMSFile(universal_constraints);

/// Interfacing with SparsePOP /// 
std::cout << "Solving with SparsePOP..." << std::endl;
std::tuple<int,int, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<double>, std::vector<domain::BoundConstraint>> fromGen(newNumVars, newNumConst, polyWidth, gndOff, gndData, observedValueById, bounds);

cp.tick("Before SparsePOP Solve");
solveWithSparsePOP(fileName, fromGen, cp);
cp.tick("After SparsePOP Solve");

int num_observations = 6;
cfg.omp_threads = std::floor(36 / num_observations);

cp.tick("Program End"); 
cp.print(); 

return 0;
}