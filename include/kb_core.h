#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <mutex>
#include <utility>
#include <algorithm>
#include <memory>

namespace kb{
//using SymID = std::uint32_t; //Do later for map-optimization
using Sym = std::string; 

struct Atom {
    Sym rel;                      
    std::vector<Sym> args;        
    std::string toString() const; 
    std::string toStringWithInput(const std::unordered_map<Sym,Sym>& freeToGround, std::unordered_map<Sym,int>& groundMap, std::vector<int>& resultVec) const; 
    // define operations for comparing Atoms
    bool operator<(const Atom &o) const noexcept;
    bool operator==(const Atom &o) const noexcept;
};

using AtomPtr = std::shared_ptr<Atom>;

using Exponent = std::uint16_t;
using MonoItem = std::pair<std::shared_ptr<Atom>, Exponent>;


struct Monomial {
    std::vector<MonoItem> items;   // kept lexicographically sorted on Atom
    void canonicalize();           // sort + merge same atom (implementation later)
    std::string toStringWithMap(const std::map<Sym, std::string>& relVarMap) const;
    std::string toStringWithInput(const std::unordered_map<Sym,Sym>& freeToGround, std::unordered_map<Sym,int>& groundMap, std::vector<int>& resultVec) const; 
    std::string toString() const;
    static std::shared_ptr<Monomial> fromAtom(const AtomPtr& a);
    static std::shared_ptr<Monomial> zeroMon();
    bool isZero() const; 
    static std::shared_ptr<Monomial> multiply(const std::shared_ptr<Monomial>& A, const std::shared_ptr<Monomial>& B);
    std::vector<AtomPtr> expandedAtoms() const;        // debug helper
    std::vector<AtomPtr> notExpandedAtoms() const; // used to generate .gms file
    bool operator<(const Monomial& o) const noexcept;  // lexicographic on items
    bool operator==(const Monomial& o) const noexcept;
};

using MonoPtr = std::shared_ptr<Monomial>;
using Coeff = double; 
using Term  = std::pair<MonoPtr, Coeff>;

struct Polynomial {
    std::vector<Term> terms;           // sorted by monomial pointer ordering
    void canonicalize();         // sort + merge same monomials (implementation later)
    static std::shared_ptr<Polynomial> fromMonomial(const MonoPtr& m);
    void addTerm(const MonoPtr& m, Coeff c);
    std::string toStringWithMap(const std::map<Sym, std::string>& relVarMap) const;
    std::string toStringWithInput(const std::unordered_map<Sym,Sym>& freeToGround, std::unordered_map<Sym,int>& groundMap, std::vector<int>& resultVec) const; 
    std::string toString() const; 
    std::string replaceString(std::string toReplace) const; 
};

enum class Cmp : std::uint8_t { EQ0, GE0 };

struct Constraint {
    Polynomial poly;
    Cmp cmp = Cmp::GE0;

    // what is this used for? 
    std::vector<std::pair<Sym,Sym>> neq;   // var‑var distinctness
    std::vector<std::string> getInputs(const std::unordered_set<Sym>& groundVariables); // Needs to be changed to support types
    void groundConstraint(std::unordered_map<Sym,int>& groundMap, const std::vector<std::string>& perm, const std::unordered_set<Sym>& groundVariables, std::string& resultString, std::vector<int>& resultVec);
    bool operator==(const Constraint& o) const noexcept;
};

}