#include "kb_core.h"
#include <iostream>
#include <stdexcept>

namespace kb {
// Atom Helpers
std::string Atom::toString() const{
    std::string out = rel + '(';
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i) out += ',';
        out += args[i];
    }
    out += ')';
    return out;
}

std::string Atom::toStringWithInput(const std::unordered_map<Sym,Sym>& freeToGround, std::unordered_map<Sym,int>& groundMap, std::vector<int>& resultVec) const{
    std::string out = rel + '(';
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i) out += ',';
        // if the argument is a free variable (exists in map made from getInputs)
        if (auto it = freeToGround.find(args[i]); it != freeToGround.end()) out += it->second; // replace it with corresponding permutation entry
        // otherwise its a groundVariable, so keep it as is
        else out += args[i];
    }
    out += ')';
    // At this point, out represents the ground Atom, so Q(x,y) is now Q(jack,jill)
    // Check if this string exists in groundMap, if not, add it
    auto it = groundMap.find(out);
    int id;
    if (it == groundMap.end()) {
        auto [iter, inserted] = groundMap.emplace(out, static_cast<int>(groundMap.size()));
        id = iter->second;           // same as groundMap.size() - 1
    } else id = it->second;             // key already present

    resultVec.push_back(id);
    return "x" + std::to_string(id);
}

bool Atom::operator<(const Atom& o) const noexcept {
    // handle zero atom (used to build zero monomial)
    if (rel.empty() && args.empty()) return true;  // empty atom is less than any other
    return std::tie(rel, args) < std::tie(o.rel, o.args);
}
bool Atom::operator==(const Atom& o) const noexcept {
    return rel == o.rel && args == o.args;
}

// Monomial Helpers
void Monomial::canonicalize() {
    std::sort(items.begin(), items.end(),
        [](const MonoItem& a, const MonoItem& b){ return *(a.first) < *(b.first); });
    std::vector<MonoItem> tmp;
    for (const auto& it : items) {
        if (!tmp.empty() && *(tmp.back().first) == *(it.first)) {
            tmp.back().second += it.second;           // same atom → add exponents
        } else {
            tmp.push_back(it);
        }
    }
    items.swap(tmp);
}

// TODO special handling on zero Monomial here
std::string Monomial::toStringWithMap(const std::map<Sym, std::string>& relVarMap) const {
    std::string out; 
    for (const auto& [ap,e] : items) {
        if (!out.empty()) out += '*';
        // Lookup atom in map
        auto it = relVarMap.find(ap->rel);
        if (it != relVarMap.end()) {
            out += it->second;  // use mapped variable name
            if (e > 1) out += "^" + std::to_string(e);
        } else {
            throw std::runtime_error("Atom is not in atom map?!?");
        }
    }
    return out; 
}
std::string Monomial::toStringWithInput(const std::unordered_map<Sym,Sym>& freeToGround, std::unordered_map<Sym,int>& groundMap, std::vector<int>& resultVec) const{
    std::string out;
    for (const auto& [ap, e] : items) {
        if (!out.empty()) out += '*';
        out += ap->toStringWithInput(freeToGround, groundMap, resultVec);
        if (e > 1) out += "^" + std::to_string(e);
    }
    return out.empty() ? "1" : out;
}

std::string Monomial::toString() const {
    std::string out;
    for (const auto& [ap, e] : items) {
        if (!out.empty()) out += '*';
        out += ap->toString();
        if (e > 1) out += "^" + std::to_string(e);
    }
    return out.empty() ? "1" : out;  // empty monomial is 1
}
std::shared_ptr<Monomial> Monomial::zeroMon(){
    auto m = std::make_shared<Monomial>();
    m->items.emplace_back(std::make_shared<Atom>(), 1); // empty atom with exponent 0
    return m;
}
bool Monomial::isZero() const {
    // a monomial is zero if it contains an empty atom with exponent 0
    return items.size() == 1 && items[0].first->rel.empty();
}
std::shared_ptr<Monomial> Monomial::fromAtom(const AtomPtr& a) {
    auto m = std::make_shared<Monomial>();
    m->items.emplace_back(a, 1);
    return m;
}

std::shared_ptr<Monomial> Monomial::multiply(const std::shared_ptr<Monomial>& A, const std::shared_ptr<Monomial>& B)
{
    // if A or B is zero monomial, then just return the other
    if (A->isZero()) return B;
    if (B->isZero()) return A; 
    // else
    auto m = std::make_shared<Monomial>();
    m->items.reserve(A->items.size() + B->items.size());
    m->items.insert(m->items.end(), A->items.begin(), A->items.end());
    m->items.insert(m->items.end(), B->items.begin(), B->items.end());
    m->canonicalize();
    return m;
}

std::vector<AtomPtr> Monomial::expandedAtoms() const {
    std::vector<AtomPtr> out;
    for (const auto& [ap, e] : items)
        for (Exponent k = 0; k < e; ++k) out.push_back(ap);
    return out;
}

std::vector<AtomPtr> Monomial::notExpandedAtoms() const {
    std::vector<AtomPtr> out;
    for (const auto& [ap, e] : items)
        out.push_back(ap);
    return out;
}


bool Monomial::operator<(const Monomial& o) const noexcept {
    // We want equality to depend on the Atoms themselves, not their pointers
    if (items.size() != o.items.size())
        return items.size() < o.items.size();

    for (std::size_t i = 0; i < items.size(); ++i) {
        const auto& [atomA, expA] = items[i];
        const auto& [atomB, expB] = o.items[i];
        if (*atomA != *atomB)         // compare atoms
            return *atomA < *atomB;
        if (expA != expB)             // and exponents
            return expA < expB;
    }
    return false; 
}
bool Monomial::operator==(const Monomial& o) const noexcept {
    if (items.size() != o.items.size()) return false;
    for (std::size_t i = 0; i < items.size(); ++i) {
        const auto& [atomA, expA] = items[i];
        const auto& [atomB, expB] = o.items[i];
        if (*atomA != *atomB) return false;
        if (expA  != expB)   return false;
    }
    return true;
}

// Polynomial Helpers //
void Polynomial::canonicalize() {
    std::sort(terms.begin(), terms.end(),
        [](const Term& a, const Term& b){ return *(a.first) < *(b.first); });
    std::vector<Term> tmp;
    for (const auto& it : terms) {
        if (!tmp.empty() && *(tmp.back().first) == *(it.first))
            tmp.back().second += it.second;        // same monomial -> add coefficients
        else
            tmp.push_back(it);
    }
    terms.swap(tmp);
}
// Create a polynomial from a single monomial 
std::shared_ptr<Polynomial> Polynomial::fromMonomial(const MonoPtr& m) {
    auto p = std::make_shared<Polynomial>();
    p->terms.emplace_back(m, 1);
    return p;
}
// something here aint working right, maybe canonicalize
void Polynomial::addTerm(const MonoPtr& m, Coeff c) {
    if (c == 0) return;
    auto cmp = [](const Term& a, const Term& b){
        return *(a.first) < *(b.first); 
        };
    auto it  = std::lower_bound(terms.begin(), terms.end(), Term{m,0}, cmp);
    if (it != terms.end() && !(*(m) < *(it->first))) {
        it->second += c;
        if (it->second == 0) terms.erase(it);
    } else {
        terms.insert(it, {m, c});
    }
    canonicalize();  // ensure polynomial is in canonical form after adding
}
 
std::string Polynomial::toStringWithMap(const std::map<Sym, std::string>& relVarMap) const {
    std::string out;
    for (const auto& [m,c] : terms) {
        if (!out.empty() && c != -1) out += " + ";
        if(m->isZero()) { // if printing zero monomial
            out += std::to_string(c); // just print coeff
            continue;
        }
        if (c == 1) {
            out += m->toStringWithMap(relVarMap);  
        } else if (c == -1) {
            out += " - " + m->toStringWithMap(relVarMap);  
        } else {
            out += std::to_string(c) + '*' + m->toStringWithMap(relVarMap); 
        } 
    }
    return out; 
}
std::string Polynomial::toStringWithInput(const std::unordered_map<Sym,Sym>& freeToGround, std::unordered_map<Sym,int>& groundMap, std::vector<int>& resultVec) const{
    std::string out; 
    for (const auto& [m,c] : terms) {
        if (!out.empty()) out += " + ";
        if(m->isZero()) { // if printing zero monomial
            out += std::to_string(c); // just print coeff
            continue;
        }
        if (c == 1) {
            out += m->toStringWithInput(freeToGround, groundMap, resultVec); 
        } else if (c == -1) {
            out += '-' + m->toStringWithInput(freeToGround, groundMap, resultVec);   
        } else {
            out += std::to_string(c) + '*' + m->toStringWithInput(freeToGround, groundMap, resultVec); 
        } 
    }
    return out;
}
std::string Polynomial::toString() const{
    std::string out; 
    for (const auto& [m,c] : terms) {
        if (!out.empty()) out += " + ";
        if(m->isZero()) { // if printing zero monomial
            //std::cout << "printing zero monomial [" << std::to_string(c) << "]" << std::endl;
            out += std::to_string(c); // just print coeff
            continue;
        }
        if (c == 1) {
            out += m->toString();  
        } else if (c == -1) {
            out += '-' + m->toString();  
        } else {
            out += std::to_string(c) + '*' + m->toString(); 
        } 
    }
    return out; 
}

std::string Polynomial::replaceString(std::string toReplace) const{
    std::string str = this->toString();
    size_t pos = 0;
    while ((pos = str.find(toReplace, pos)) != std::string::npos) {
        str.replace(pos, toReplace.length(), "G");
        pos += 1;
    }
    return str;
}

static bool termVecEqual(const std::vector<Term>& a, const std::vector<Term>& b){
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].second != b[i].second) return false; // coeff
        if (a[i].first->toString() != b[i].first->toString()) // monomial
            return false;
    }
    return true;
}

bool Constraint::operator==(const Constraint& o) const noexcept
{
    return cmp == o.cmp
        && neq == o.neq
        && termVecEqual(poly.terms, o.poly.terms);
}

// return all free variables in a constraint
std::vector<std::string> Constraint::getInputs(const std::unordered_set<Sym>& groundVariables){
    std::unordered_set<std::string> inputSet; 
    for (const auto& term : poly.terms) { // each term is a monomial
        for (const auto& monoItem : term.first->items) {
            for (const auto& arg : monoItem.first->args) {
                // if arg is not in groundVariables, add it to inputs
                if (groundVariables.find(arg) == groundVariables.end()) {
                    inputSet.insert(arg);
                }
            }
        }
    }
    std::vector<std::string> inputs;
    for (auto& var : inputSet) inputs.push_back(var);
    return inputs;
}

std::vector<std::pair<SymbolType, std::string>> Constraint::getOrderedTypedInputs() const {
    std::vector<std::pair<SymbolType, std::string>> result; 
    std::unordered_set<std::string> seen; 

    for (const auto& term : poly.terms) { 
        for (const auto& monoItem : term.first->items) {
            for (const auto& arg : monoItem.first->args) {
                if (seen.find(arg) != seen.end()) continue; 

                SymbolType type; 
                if (arg.find("gene") == 0) {
                    type = SymbolType::GENE;
                } else if (arg.find("enzyme") == 0) {
                    type = SymbolType::ENZYME;
                } else if (arg.find("reaction") == 0) {
                    type = SymbolType::REACTION;
                } else if (arg.find("compound") == 0) {
                    type = SymbolType::COMPOUND;
                } else {
                    throw std::runtime_error("Unknown symbol type for argument: " + arg);
                }

                result.push_back({type, arg}); 
                seen.insert(arg); 
            }
        }
    }
    return result; 
}

std::vector<int> Constraint::groundToAtomIDs(const std::unordered_map<Sym,std::string>& substitution, std::unordered_map<size_t,int>& groundMap) const {
    std::vector<int> atomIDs; 
    for (const auto& term : poly.terms) {
        // skip over zero terms -> constants in constraints like 1 in "friends(x,y)-1>0"
        if (term.first->isZero()) continue;

        for (const auto& monoItem : term.first->items) {
            const auto& atom = monoItem.first; 
            std::string groundAtomStr = atom->rel + '(';

            for (std::size_t i = 0; i < atom->args.size(); i++) {
                const auto& arg = atom->args[i];

                auto it = substitution.find(arg); 
                if (it == substitution.end()) {
                    throw std::runtime_error("Substitution missing argument: " + arg); 
                }
                groundAtomStr += it->second;
                // add comma between arguments but not after last one
                if (i < atom->args.size() -1) {
                    groundAtomStr += ","; 
                }
            }
            groundAtomStr += ")"; 
            // Hash resulting string and look it up in groundMap
            size_t hashValue = std::hash<std::string>{}(groundAtomStr);
            auto mapIt = groundMap.find(hashValue); 
            int atomID; 

            if (mapIt == groundMap.end()) {
                // ground atom hasn't been seen before
                // assign it to next available ID
                atomID = static_cast<int>(groundMap.size());
                groundMap[hashValue] = atomID; 
            } else {
                atomID = mapIt->second; 
            }
            // if (atomID%1000 == 0 ) { std::cout << "AtomID: " << atomID << std::endl; }
            //std::cout << "Ground Atom: " << groundAtomStr << " -> AtomID: " << atomID << std::endl;
            atomIDs.push_back(atomID);
        }
    }
    return atomIDs; 
}

void Constraint::groundConstraint(std::unordered_map<Sym,int>& groundMap, const std::vector<std::string>& perm, const std::unordered_set<Sym>& groundVariables, std::string& resultString, std::vector<int>& resultVec){
    std::vector<std::string> inputs = getInputs(groundVariables);
    if (inputs.size() != perm.size()){ 
        std::cerr << "Trying to ground constraint that takes " << inputs.size() << " variables with a permutation of " << perm.size() << " variables!" << std::endl; 
        return; 
    }

    // freeToGround holds a map from free variables in constraint to ground variables in permutation
    std::unordered_map<Sym, Sym> freeToGround; 
    for (std::size_t i = 0; i < inputs.size(); i++){ // populate freeToGround
        freeToGround[inputs[i]] = perm[i];
    }

    for(auto& cond : neq){ // Check that is satisfies conditions like x!=y
        Sym first, second; 
        if (freeToGround.find(cond.first) != freeToGround.end()) first = freeToGround[cond.first]; // free variable
        else first = cond.first; // if ground variable
        if (freeToGround.find(cond.second) != freeToGround.end()) second = freeToGround[cond.second]; // free variable
        else second = cond.second; // if ground variable
        if (first == second) return; 
    }

    // Print out new constraint
    resultString = poly.toStringWithInput(freeToGround, groundMap, resultVec);
}

}