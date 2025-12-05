#ifndef DOMAIN_H
#define DOMAIN_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "kb_core.h"

namespace domain {

using kb::SymbolType;
using kb::Sym;
using kb::Coeff;
using kb::AtomPtr; 
using kb::Atom;
using kb::MonoPtr;
using kb::Monomial;
using kb::Polynomial;
using kb::Cmp;
using kb::Constraint;

// Types of ground symbols in our domain
//enum class SymbolType : std::uint8_t {GENE, ENZYME, REACTION, COMPOUND};

// Global map of predicate names to their argument type signatures.
inline std::unordered_map<std::string, std::vector<SymbolType>> PREDICATE_SIGNATURES;

struct GroundNames { // stores TYPED ground names
    std::unordered_set<std::string> genes;
    std::unordered_set<std::string> enzymes;
    std::unordered_set<std::string> reactions;
    std::unordered_set<std::string> compounds;
};

// Initialize the global predicate signatures map
void initializePredicateSignatures();

// For parsing ProbLog Files, ex: 
// 0.183::function(g614,ec_4_13).
// ortholog(g614,g616).             (this is implicitly 1.0::ortholog(g614,g616).)
class ProbLogParser {
public:
    ProbLogParser(GroundNames& gn);
    std::vector<kb::Constraint> parseFile(const std::string& filename); 
private:
    GroundNames& groundNames;

    // parse one line -> get one constraint
    kb::Constraint parseLine(const std::string& line);
    // extract prob and atom from a single line
    std::pair<double, std::string> extractProbAndAtom(const std::string& line);
    // extract probability
    kb::Atom parseAtom(const std::string& atomStr);
    // Build simple constraint: atom - prob = 0
    kb::Constraint buildConstraint(const kb::Atom& atom, double prob);

    // Helpers
    std::vector<std::string> splitArgs(const std::string& argsStr); // split by comma
    std::string trim(const std::string& s); // trim whitespace 
};

Constraint parseConstraint(const std::string &text);

void generateGrounding(const std::vector<kb::Constraint>& constraints, const std::vector<std::vector<std::string>>& typedGroundNames, std::unordered_map<size_t,int>& groundMap, std::vector<std::vector<std::vector<int>>>& resultVec);

void createGroundingRepresentation(const std::vector<std::vector<std::vector<int>>>& finalResults, std::vector<int>& polyWidth, std::vector<int>& gndOff, std::vector<int>& gndData);
// std::map<std::string, std::string> relVarMap(const std::vector<kb::Constraint>& constraints);
std::string writeGMSFile(const std::vector<kb::Constraint>& constraints); 

} 
#endif