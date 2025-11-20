#ifndef DOMAIN_H
#define DOMAIN_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include "kb_core.h"

namespace domain {

// Types of ground symbols in our domain
enum class SymbolType : std::uint8_t {GENE, ENZYME, REACTION, COMPOUND};

// Global map of predicate names to their argument type signatures.
inline std::unordered_map<std::string, std::vector<SymbolType>> PREDICATE_SIGNATURES;

struct GroundNames {
    std::unordered_set<std::string> genes;
    std::unordered_set<std::string> enzymes;
    std::unordered_set<std::string> reactions;
    std::unordered_set<std::string> compounds;
};

// Initialize the global predicate signatures map
void initializePredicateSignatures();

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

} 
#endif