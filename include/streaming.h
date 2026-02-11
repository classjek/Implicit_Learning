#ifndef _STREAMING_H_
#define _STREAMING_H_

#include <vector>
#include <unordered_map>
#include <cstdio>
#include <functional>
#include <string>

// Forward declarations
class spvec_array;
class poly_info;

// MonomialKey: Represents a monomial for hashing (sparse: var_index -> exponent)
struct MonomialKey {
    std::vector<std::pair<int, int>> terms;  // (variable_index, exponent) pairs, sorted by index
    MonomialKey() = default;
    
    // Construct from spvec_array at a specific monomial position
    MonomialKey(const spvec_array& sups, int monomial_idx);
    bool operator==(const MonomialKey& other) const {
        return terms == other.terms;
    }
};

// Hash function for MonomialKey
struct MonomialKeyHash {
    size_t operator()(const MonomialKey& key) const {
        size_t h = 0;
        for (const auto& [var, exp] : key.terms) {
            // Combine hashes: var and exponent
            h ^= std::hash<int>()(var) + 0x9e3779b9 + (h << 6) + (h >> 2);
            h ^= std::hash<int>()(exp) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

// StreamingContext: Shared state for both passes
struct StreamingContext {
    // Monomial -> variable number mapping (built in pass 1, used in pass 2)
    std::unordered_map<MonomialKey, int, MonomialKeyHash> monomial_to_var;
    
    // Block structure info
    std::vector<int> block_struct; // Size of each block (positive = matrix, negative = diagonal)
    std::vector<int> block_offsets; // Starting entry index for each block
    
    int mDim = 0; // Total number of SDP variables (unique monomials)
    int nBlocks = 0; // Total number of blocks
    
    // Objective coefficients (indexed by variable number)
    std::vector<double> obj_coef;
    // used in pass 2
    FILE* output_file = nullptr;
    // Pass indicator
    bool is_counting_pass = true;
    // Current block being processed
    int current_block = 0;
    // Entry count for current block
    int current_block_entries = 0;
    int total_entries = 0; // total count of matrix entries written 

    const std::vector<int>* binvec_ptr = nullptr;
    const std::vector<int>* Sqvec_ptr = nullptr;
    
    
    // Register a monomial, returns its variable number (1-indexed for SDPA format)
    int register_monomial(const MonomialKey& key);
    
    // Look up variable number for a monomial (pass 2)
    int get_var_number(const MonomialKey& key) const;
    
    // Start a new block
    void start_block(int block_size);
    
    // Write an SDP entry (pass 2 only)
    void write_entry(int var_num, int block, int row, int col, double coef);
    
    // Finalize pass 1 (prepare for pass 2)
    void finalize_counting();
    
    // Write SDPA header (beginning of pass 2)
    void write_header(const std::string& filename);
    
    // Write SDPA footer / finalize file
    void finalize_file();
};


// streaming functions
void convert_obj_stream(class poly_info& polyinfo, StreamingContext& ctx);
void convert_eq_stream(class poly_info& polyinfo, class spvec_array& bassinfo, StreamingContext& ctx);
void convert_ineq_a_ba1_stream(class poly_info& polyinfo, class spvec_array& bassinfo, StreamingContext& ctx);
void convert_ineq_a_ba2_stream(class poly_info& polyinfo, class spvec_array& bassinfo, StreamingContext& ctx);
void convert_sdp_stream(class poly_info& polyinfo, class spvec_array& bassinfo, StreamingContext& ctx);
void convert_ba1mmt_stream(class spvec_array& bassinfo, StreamingContext& ctx);
void convert_ba2mmt_stream(class spvec_array& bassinfo, StreamingContext& ctx);

// replaces get_psdp + write_sdpa
void stream_psdp_to_file(
    int mdim,
    int msize,
    std::vector<class poly_info>& polyinfo,
    std::vector<class spvec_array>& bassinfo,
    const std::string& sdpafile, 
    const std::vector<int>& binvec = {}, 
    const std::vector<int>& Sqvec = {}
);

void test_streaming_basics(); 

#endif // _STREAMING_H_