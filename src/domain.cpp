#include "domain.h"
#include "kb_core.h"
#include <fstream>
#include <stdexcept>

namespace domain {

void initializePredicateSignatures() {
    // Custom to our domain
    PREDICATE_SIGNATURES["function"] = {SymbolType::GENE, SymbolType::ENZYME};
    PREDICATE_SIGNATURES["ortholog"] = {SymbolType::GENE, SymbolType::GENE};
    PREDICATE_SIGNATURES["reaction_enzyme"] = {SymbolType::REACTION, SymbolType::ENZYME};
    PREDICATE_SIGNATURES["reaction_compound_reaction"] = {SymbolType::REACTION, SymbolType::COMPOUND, SymbolType::REACTION};
    PREDICATE_SIGNATURES["accept_compound"] = {SymbolType::COMPOUND};
}

// Add Problog Parser 
ProbLogParser::ProbLogParser(GroundNames& gn) : groundNames(gn) {}

std::vector<kb::Constraint>ProbLogParser::parseFile(const std::string& filename) {
    std::vector<kb::Constraint> constraints;
    std::ifstream file(filename); 

    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::string line; 
    while (std::getline(file, line)) {
        line = trim(line); 

        if (line.empty() || line[0] == '%') continue; // skip empty lines and comments
        kb::Constraint constraint = parseLine(line); 
        constraints.push_back(constraint); 
    }
    return constraints;
}

kb::Constraint ProbLogParser::parseLine(const std::string& line) {
    // Extract probability and atom string
    auto [prob, atomStr] = extractProbAndAtom(line);
    
    // Parse the atom (and register ground names)
    kb::Atom atom = parseAtom(atomStr);
    
    // Build constraint: atom - prob = 0
    return buildConstraint(atom, prob);
}

std::pair<double, std::string> ProbLogParser::extractProbAndAtom(const std::string& line) {
    double prob = 1.0; 
    std::string atomStr; 

    size_t colonPos = line.find("::"); 
    if (colonPos != std::string::npos) {
        // Probability specified, ex: 0.183::function(g614,ec_4_13).
        std::string probStr = line.substr(0, colonPos); 
        prob = std::stod(probStr); 
        atomStr = line.substr(colonPos + 2); 
    } else { 
        // no associated prob, default to 1.0. ex: accept_compound(c16875). 
        atomStr = line; 
    }

    // TODO: is this double checking needed? 
    //remove trailing '.'
    atomStr = trim(atomStr);
    if (!atomStr.empty() && atomStr.back() == '.') {
        atomStr.pop_back(); 
    }
    return {prob, atomStr}; 
}

kb::Atom ProbLogParser::parseAtom(const std::string& atomStr) {
    size_t parenPos = atomStr.find('('); 
    if (parenPos == std::string::npos) {
        throw std::runtime_error("Invalid atom format: " + atomStr); 
    }

    // extract ground name
    std::string predicate = atomStr.substr(0, parenPos); 
    size_t closeParenPos = atomStr.find(')', parenPos); 
    if (closeParenPos == std::string::npos) {
        throw std::runtime_error("Invalid atom format (no closing parenthesis): " + atomStr);
    }

    // extract arguments
    size_t argStart = parenPos + 1; 
    size_t argLen = closeParenPos - argStart; 
    std::string argsStr = atomStr.substr(argStart, argLen); 

    // split arguments by comma 
    std::vector<std::string> args = splitArgs(argsStr); 
    // Look up predicate signature
    auto it = PREDICATE_SIGNATURES.find(predicate); 
    if (it == PREDICATE_SIGNATURES.end()) {
        throw std::runtime_error("Unknown predicate: " + predicate); 
    }

    const std::vector<SymbolType>& signature = it->second;
    
    //validate arg count
    if (args.size() != signature.size()) {
        throw std::runtime_error("Argument count mismatch for predicate " + predicate); 
    }

    // register (typed) ground names
    for (size_t i = 0; i < args.size(); i++) {
        std::string arg = trim(args[i]); 
        SymbolType type = signature[i]; 

        switch(type) { 
            case SymbolType::GENE:
                groundNames.genes.insert(arg); 
                break; 
            case SymbolType::ENZYME:
                groundNames.enzymes.insert(arg); 
                break; 
            case SymbolType::REACTION:
                groundNames.reactions.insert(arg); 
                break; 
            case SymbolType::COMPOUND:
                groundNames.compounds.insert(arg); 
                break;
        }
        args[i] = arg;
    }

    // Build and return atom
    kb::Atom atom;
    atom.rel = predicate;
    atom.args = args;
    
    return atom;
}

kb::Constraint ProbLogParser::buildConstraint(const kb::Atom& atom, double prob) {
    kb::Constraint constraint; 

    // Add Atom with coefficient 1.0
    auto atomPtr = std::make_shared<kb::Atom>(atom);
    auto monoPtr = kb::Monomial::fromAtom(atomPtr); 
    constraint.poly.addTerm(monoPtr, 1.0); 

    // add probability 
    auto zeroMono = kb::Monomial::zeroMon(); 
    constraint.poly.addTerm(zeroMono, -prob);

    constraint.poly.canonicalize(); 
    // set comparison to equality
    constraint.cmp = kb::Cmp::EQ0;

    return constraint; 
}

std::vector<std::string> ProbLogParser::splitArgs(const std::string& argsStr) {
    std::vector<std::string> result;
    std::string current;
    
    for (char c : argsStr) {
        if (c == ',') {
            result.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    // Add last argument
    if (!current.empty()) {
        result.push_back(current);
    }
    
    return result;
}

std::string ProbLogParser::trim(const std::string& s) {
    size_t start = 0;
    size_t end = s.length();
    // Find first non-whitespace
    while (start < end && std::isspace(s[start])) {
        ++start;
    }
    //find last non-whitespace
    while (end > start && std::isspace(s[end - 1])) {
        --end;
    }
    
    return s.substr(start, end - start);
}

}