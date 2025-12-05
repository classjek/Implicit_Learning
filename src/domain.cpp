#include "domain.h"

#include <fstream>
#include <functional> 
#include <iostream>
#include <stdexcept>
#include <unordered_map>

#include "kb_core.h"

namespace domain {

void initializePredicateSignatures() { // custom to our domain
    PREDICATE_SIGNATURES["function"] = {SymbolType::GENE, SymbolType::ENZYME};
    PREDICATE_SIGNATURES["ortholog"] = {SymbolType::GENE, SymbolType::GENE};
    PREDICATE_SIGNATURES["reaction_enzyme"] = {SymbolType::REACTION, SymbolType::ENZYME};
    PREDICATE_SIGNATURES["reaction_compound_reaction"] = {SymbolType::REACTION, SymbolType::COMPOUND, SymbolType::REACTION};
    PREDICATE_SIGNATURES["accept_compound"] = {SymbolType::COMPOUND};
    PREDICATE_SIGNATURES["reaction"] = {SymbolType::REACTION, SymbolType::COMPOUND, SymbolType::REACTION}; 
    PREDICATE_SIGNATURES["enzyme_reaction_path"] = {SymbolType::GENE, SymbolType::ENZYME, SymbolType::REACTION, SymbolType::REACTION, SymbolType::ENZYME, SymbolType::GENE};
    PREDICATE_SIGNATURES["ortholog_support"] = {SymbolType::GENE, SymbolType::GENE, SymbolType::ENZYME}; 
}

#pragma region ProbLogParser
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

#pragma endregion // end ProbLogParser

#pragma region ConstraintParser 
// The ConstraintParser is only exposed via the domain.h interface function parseConstraint
// thus we implement it entirely here in domain.cpp
enum class Tok { IDENT, NUMBER, PLUS, MINUS, STAR, LP, RP, COMMA, GE, EQ, NEQ, COLON, END };
struct Token { Tok kind; std::string text; };

std::unordered_map<std::string, AtomPtr> atomPool;
std::unordered_map<std::string, MonoPtr> monomialPool;

// canonical-string builder (same as your makeKey)
static std::string atomKey(const std::string& rel, const std::vector<Sym>& args) {
    std::string k = rel;
    k.push_back('(');
    for (std::size_t i=0;i<args.size();++i) {
        k += args[i];
        if (i+1<args.size()) k.push_back(',');
    }
    k.push_back(')');
    return k;
}

AtomPtr internAtom(const std::string& rel, const std::vector<Sym>& args) {
    std::string key = atomKey(rel, args);
    auto [it, inserted] = atomPool.try_emplace(key, nullptr);
    if (inserted) {
        auto a = std::make_shared<Atom>();
        a->rel  = rel;
        a->args = args;
        it->second = a;
    }
    return it->second;
}

class Lexer {
public:
    explicit Lexer(const std::string &s) : src(s), p(src.c_str()) { next(); }

    const Token &peek() const { return cur; }
    Token pop()  { Token t = cur; next(); return t; }
    Token peekNext() const {
        //const char *saveP = p;
        Token saveCur     = cur;
        // create a throw-away lexer just to advance once
        Lexer tmp(*this);
        tmp.pop();                 // skip current token
        Token nxt = tmp.peek();
        return nxt;
    }

private:
    void next() {
        while (std::isspace(*p)) ++p;
        if (*p == '\0') { cur = {Tok::END, ""}; return; }
        char c = *p;
        if (std::isalpha(c) || c=='_') {
            const char *start = p;
            while (std::isalnum(*p) || *p=='_') ++p;
            cur = {Tok::IDENT, std::string(start, p)}; return;
        }
        if (std::isdigit(c) || (c == '.' && std::isdigit(*(p + 1)))) {
            const char* start = p;
            bool seenDot = false;
            while (std::isdigit(*p) || (*p == '.' && !seenDot)) {
                if (*p == '.') seenDot = true;
                ++p;
            }
            cur = {Tok::NUMBER, std::string(start, p)};
            return;
        }
        ++p; // consume single char tokens
        switch (c) {
            case '+': cur = {Tok::PLUS, "+"}; return;
            case '-': cur = {Tok::MINUS,"-"}; return;
            case '*': cur = {Tok::STAR, "*"}; return;
            case '(': cur = {Tok::LP,   "("}; return;
            case ')': cur = {Tok::RP,   ")"}; return;
            case ',': cur = {Tok::COMMA,","}; return;
            case ':': cur = {Tok::COLON, ":"}; return;
            case '>': if (*p=='='){ ++p; cur={Tok::GE, ">="}; return; }; break;
            case '!': if (*p=='='){ ++p; cur = {Tok::NEQ, "!="}; return; } break;
            default:  break;
        }
        if (c=='=') { cur={Tok::EQ, "="}; return; }
        throw std::runtime_error("Unexpected char in input");
    }
    const std::string &src; const char *p; Token cur;
};


struct Parser {
    explicit Parser(const std::string &s) : lex(s) {}
    Constraint parse();
private:
    // grammar helpers
    Polynomial parseSum();              // sum  ::= product ((+|-) product)*
    std::pair<bool,Coeff>  parseCoefficient();      // coefficient ::= number? '*'
    MonoPtr    parseProduct();          // product ::= factor (factor)*
    MonoPtr    parseFactor();           // factor ::= atom | number
    AtomPtr    parseAtom();

    // utility
    bool accept(Tok k){ if(lex.peek().kind==k){ lex.pop(); return true;} return false; }
    void expect(Tok k,const char*msg){ if(!accept(k)) throw std::runtime_error(msg);}  
    std::size_t varIndex(const std::string& name, std::vector<std::string>& vars);

    Lexer lex;
};

std::size_t varIndex(const std::string& name, std::vector<std::string>& vars) {
    auto it = std::find(vars.begin(), vars.end(), name);
    if (it != vars.end()) return std::distance(vars.begin(), it);
    vars.push_back(name);
    return vars.size()-1;
}

AtomPtr Parser::parseAtom() {
    // Extract id and args, i.e. Q and x,b from Q(x,b)
    Token id = lex.pop();              
    expect(Tok::LP, "Expected '('");
    std::vector<Sym> args;
    if (lex.peek().kind != Tok::RP) {
        do {
            if (lex.peek().kind != Tok::IDENT)
                throw std::runtime_error("Expected identifier in arg list");
            args.push_back(lex.pop().text);
        } while (accept(Tok::COMMA));
    }
    expect(Tok::RP, "Expected ')'");

    // Add atom to canonical atom store if not already there
    // return pointer to atom in store
    return internAtom(id.text, args);
}

MonoPtr Parser::parseFactor() {
    if (lex.peek().kind == Tok::IDENT) {
        return Monomial::fromAtom(parseAtom());
    }
    if (lex.peek().kind == Tok::NUMBER) {
        // treat numeric constant n as n * 1, handled in Polynomial layer
        int val = std::stoi(lex.pop().text);
        auto one = std::make_shared<Monomial>(); // empty product == 1
        auto monoPtr = one;
        Polynomial dummy;
        dummy.addTerm(monoPtr, val);            // store constant in dummy
        // pull the monomial pointer back out (same obj) for caller to use
        return monoPtr;                         // coeff handled in parseProduct
    }
    // throw std::runtime_error("Unexpected token in factor");
    throw std::runtime_error(
        "Unexpected token in factor: kind=" +
        std::to_string(static_cast<int>(lex.peek().kind)) +
        " text='" + lex.peek().text + "'");
}

std::pair<bool, Coeff> Parser::parseCoefficient() {
    if (lex.peek().kind != Tok::NUMBER)
        return {true, static_cast<Coeff>(1.0)};      // Coeff is now double

    std::string numStr = lex.pop().text;
    Coeff value = static_cast<Coeff>(std::stod(numStr));   // handles 42, 0.3, .25 …

    if (lex.peek().kind != Tok::STAR)
        return {false, value};

    accept(Tok::STAR);
    return {true, value};
}


MonoPtr Parser::parseProduct() {
    auto m = parseFactor();
    // Implicit multiplication: IDENT IDENT … or with '*'
    while (lex.peek().kind == Tok::IDENT || lex.peek().kind == Tok::LP || lex.peek().kind == Tok::STAR) {
        if (accept(Tok::STAR)) continue;        // consume explicit '*'
        auto rhs = parseFactor();
        m = Monomial::multiply(m, rhs);
    }
    return m;
}

Polynomial Parser::parseSum() {
    Polynomial P;
    bool neg = false; // handle optional leading sign (+/-)
    if (accept(Tok::PLUS) || (neg = accept(Tok::MINUS))) {}
    
    // parse first object and add its monomial to polynomial
    auto coef = parseCoefficient();
    if (coef.first){ // there is following term after coeff
        auto firstMono = parseProduct(); // process term
        P.addTerm(firstMono, neg ? -1 * coef.second : 1 * coef.second);
    } else { // no following term 
        // add zeroMonomial to represent constant 
        P.addTerm(Monomial::zeroMon(), neg ? -1 * coef.second : 1 * coef.second);
    }

    // add remaining terms in the sum
    while (lex.peek().kind == Tok::PLUS || lex.peek().kind == Tok::MINUS) {
        neg = accept(Tok::MINUS); 
        if (!neg) accept(Tok::PLUS);
        auto coef = parseCoefficient(); 
        if (coef.first){ // there is term after coeff
            auto m = parseProduct(); // process term 
            P.addTerm(m, neg ? -1 * coef.second : 1 * coef.second);
        } else { // no following term
            P.addTerm(Monomial::zeroMon(), neg ? -1 * coef.second : 1 * coef.second);
        }
    }
    return P;
}

Constraint Parser::parse() {
    Constraint C;
    // handle distinctness guard like x != y : 
    std::vector<std::string> vars;               // keeps order of seen variables

    // Handle distinctness guard: x != y
    if (lex.peek().kind == Tok::IDENT && lex.peekNext().kind == Tok::NEQ) {
        do {
            Token a = lex.pop();                // IDENT
            expect(Tok::NEQ, "need '!=' in guard");
            Token b = lex.pop();                // IDENT

            C.neq.emplace_back(a.text, b.text);

        } while (accept(Tok::COMMA));
        expect(Tok::COLON, "missing ':' after guard");
    }

    // handles the left hand side of the constraint, generating its polynomial representation
    Polynomial lhs = parseSum();

    Tok compTok = lex.peek().kind;
    if (compTok == Tok::GE || compTok == Tok::EQ) lex.pop();
    else throw std::runtime_error("Expected '>=' or '='");

    Polynomial rhs = parseSum();

    // Move rhs to lhs just in case 
    for (auto [m,c] : rhs.terms) lhs.addTerm(m, -c);
    C.poly = std::move(lhs);
    C.cmp  = (compTok == Tok::EQ ? Cmp::EQ0 : Cmp::GE0);

    if (lex.peek().kind != Tok::END)
        throw std::runtime_error("Unexpected trailing tokens");

    return C;
}

// Public API
Constraint parseConstraint(const std::string &text) {
    Parser p(text);
    return p.parse();
}

#pragma endregion // end ConstraintParser

#pragma region Grounding

// void groundConstraint(const std::vector<std::pair<SymbolType, std::string>>& grounding, std::unordered_map<Sym,int>& groundMap, std::vector<std::vector<std::vector<int>>>& resultVec) {
//     // placeholder implementation
//     //std::cout << "Called groundConstraint with grounding size: " << grounding.size() << std::endl;
//     for (const auto& [symType, name] : grounding) {
//         // std::cout << static_cast<int>(symType) << ": " << name << ", ";
//         std::cout << name << ", ";
//     }
//     std::cout << std::endl;
//     // Actual implementation would map the grounding to indices and store in resultVec
// }
void groundConstraint(const kb::Constraint& constraint, const std::vector<std::pair<SymbolType, std::string>>& orderedTypedInputs, 
     const std::vector<std::pair<SymbolType, std::string>>& grounding, std::unordered_map<size_t, int>& groundMap, std::vector<std::vector<int>>& constraintGroundings) {

        std::unordered_map<SymbolType, int> typeCounter; 
        std::unordered_map<std::string, std::string> substitution; 

        for (const auto& [type, varName] : orderedTypedInputs) {
            // skip if we've already mapped this variable
            if (substitution.find(varName) != substitution.end()) {
                continue; 
            }

            int targetIdx = typeCounter[type]; 
            int currentIdx = 0; 

            for (const auto& [gType, gName] : grounding) {
                if (gType == type) {
                    if (currentIdx == targetIdx) {
                        substitution[varName] = gName; 
                        break;
                    }
                    currentIdx++; 
                }
            }
            typeCounter[type]++;
        }
        std::vector<int> atomIDs = constraint.groundToAtomIDs(substitution, groundMap);
        constraintGroundings.push_back(atomIDs);
}

// Later on, this can be made parallel per constraint
// but DFS is the main bottleneck here, so that will only help if we have many constraints
void generateGrounding(const std::vector<kb::Constraint>& constraints, const std::vector<std::vector<std::string>>& typedGroundNames, std::unordered_map<size_t,int>& groundMap, std::vector<std::vector<std::vector<int>>>& resultVec) {
    std::cout << "Called generateGrounding on all Constraints" << std::endl;
    resultVec.resize(constraints.size());

    // for all constraints
    for (size_t i = 0; i < constraints.size(); i++) {
        // get input information for this constraint
        std::vector<std::pair<SymbolType, std::string>> orderedTypedInputs = constraints[i].getOrderedTypedInputs(); 
        std::vector<std::pair<SymbolType, int>> typeSequence;
        std::unordered_map<SymbolType,int> countMap; 

        for (const auto& [symType, name] : orderedTypedInputs) {
            countMap[symType]++;
        }
        for (const auto& [symType, count] : countMap) {
            typeSequence.push_back({symType, count}); 
        } 

        // grounding that we'll build up 
        std::vector<std::pair<SymbolType, std::string>> grounding; 

        // nested DFS function to generate all type-aware groundings
        std::function<void(int, int)> dfs = [&](int typeIdx, int countRemaining) -> void {
        //  Base case 1: finished with all types - have full grounding 
        if (typeIdx >= static_cast<int>(typeSequence.size())) {
            //groundConstraint(grounding, groundMap, resultVec);
            groundConstraint(constraints[i], orderedTypedInputs, grounding, groundMap, resultVec[i]);  
            return; 
        }
        // Base case 2: finished with current type - move to next type
        if (countRemaining == 0) {
            dfs(typeIdx + 1, typeSequence[typeIdx + 1].second);
            return;
        }
        // recursive case: pick a name of the current type
        SymbolType currentType = typeSequence[typeIdx].first; 

        // Loop through all available names for this type
        for (const std::string& name : typedGroundNames[static_cast<int>(currentType)]) {
            // add to grounding 
            grounding.push_back({currentType, name}); 
            // recurse with one less to pick of this type
            dfs(typeIdx, countRemaining - 1);
            grounding.pop_back(); 
            }
        };


        if (!typeSequence.empty()) {
            dfs(0, typeSequence[0].second);
        }
    }
}

void createGroundingRepresentation(const std::vector<std::vector<std::vector<int>>>& finalResults, std::vector<int>& polyWidth, std::vector<int>& gndOff, std::vector<int>& gndData) {
    // make up for dummy objective function: first constraint, takes no arguments
    polyWidth.push_back(0);
    gndOff.push_back(0);

    for (const auto& constraint : finalResults) {
        // add number of arguments taken for given constraint
        polyWidth.push_back(constraint.empty() ? 0 : static_cast<int>(constraint[0].size()));
        // offset for this constraint
        gndOff.push_back(static_cast<int>(gndData.size()));
        // add all groudings for this constraint
        for (const auto& grounding : constraint) {
            for (int atomID : grounding) {
                gndData.push_back(atomID);
            }
        }
    }
}


#pragma endregion
}