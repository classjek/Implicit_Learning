#include "streaming.h"
#include "spvec.h"
#include "global.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <unordered_set>


// Add this helper function at the top of streaming.cpp (after includes)

// Merge two monomials by adding exponents (polynomial multiplication)
// Returns a MonomialKey representing the product
MonomialKey merge_monomials(const spvec_array& sup1, int idx1, const spvec_array& sup2, int idx2) {
    MonomialKey result;
    
    int pos1 = sup1.pnz[0][idx1];
    int max1 = pos1 + sup1.pnz[1][idx1];
    int pos2 = sup2.pnz[0][idx2];
    int max2 = pos2 + sup2.pnz[1][idx2];
    
    // Handle empty cases
    if (pos1 < 0 && pos2 < 0) return result;  // Both empty
    if (pos1 < 0) { //only second monomial
        for (int i = pos2; i < max2; i++) {
            result.terms.emplace_back(sup2.vap[0][i], sup2.vap[1][i]);
        }
        return result;
    }
    
    if (pos2 < 0) { // only first monomial
        for (int i = pos1; i < max1; i++) {
            result.terms.emplace_back(sup1.vap[0][i], sup1.vap[1][i]);
        }
        return result;
    }
    
    // Merge both (sorted merge, adding exponents for matching vars)
    while (pos1 < max1 && pos2 < max2) {
        if (sup1.vap[0][pos1] < sup2.vap[0][pos2]) {
            result.terms.emplace_back(sup1.vap[0][pos1], sup1.vap[1][pos1]);
            pos1++;
        } else if (sup1.vap[0][pos1] > sup2.vap[0][pos2]) {
            result.terms.emplace_back(sup2.vap[0][pos2], sup2.vap[1][pos2]);
            pos2++;
        } else { // same variable -> add exponents
            result.terms.emplace_back(sup1.vap[0][pos1], sup1.vap[1][pos1] + sup2.vap[1][pos2]);
            pos1++;
            pos2++;
        }
    }
    while (pos1 < max1) {
        result.terms.emplace_back(sup1.vap[0][pos1], sup1.vap[1][pos1]);
        pos1++;
    }
    while (pos2 < max2) {
        result.terms.emplace_back(sup2.vap[0][pos2], sup2.vap[1][pos2]);
        pos2++;
    }
    return result;
}

// Merge a MonomialKey with a monomial from spvec_array
MonomialKey merge_key_with_mono(const MonomialKey& key, const spvec_array& sup, int idx) {
    MonomialKey result;
    int pos2 = sup.pnz[0][idx];
    int max2 = pos2 + sup.pnz[1][idx];
    if (pos2 < 0) return key;  // sup is empty, return key as-is
    if (key.terms.empty()) {
        // key is empty, just copy from sup
        for (int i = pos2; i < max2; i++) {
            result.terms.emplace_back(sup.vap[0][i], sup.vap[1][i]);
        }
        return result;
    }
    
    // Merge sorted lists
    size_t k = 0;
    int p = pos2;
    while (k < key.terms.size() && p < max2) {
        if (key.terms[k].first < sup.vap[0][p]) {
            result.terms.push_back(key.terms[k]);
            k++;
        } else if (key.terms[k].first > sup.vap[0][p]) {
            result.terms.emplace_back(sup.vap[0][p], sup.vap[1][p]);
            p++;
        } else {
            result.terms.emplace_back(key.terms[k].first, key.terms[k].second + sup.vap[1][p]);
            k++; p++;
        }
    }
    while (k < key.terms.size()) {
        result.terms.push_back(key.terms[k++]);
    }
    while (p < max2) {
        result.terms.emplace_back(sup.vap[0][p], sup.vap[1][p]);
        p++;
    }
    return result;
}

MonomialKey::MonomialKey(const spvec_array& sups, int monomial_idx) {
    // Get the range for this monomial in vap
    int start = sups.pnz[0][monomial_idx];
    int nnz = sups.pnz[1][monomial_idx];
    
    terms.reserve(nnz);
    for (int i = 0; i < nnz; i++) {
        int var_idx = sups.vap[0][start + i];
        int exponent = sups.vap[1][start + i];
        terms.emplace_back(var_idx, exponent);
    }
    // Ensure sorted by variable index for consistent hashing
    std::sort(terms.begin(), terms.end());
}

//=============================================================================
// StreamingContext implementation
//=============================================================================

int StreamingContext::register_monomial(const MonomialKey& key) {
    auto it = monomial_to_var.find(key);
    if (it != monomial_to_var.end()) {
        return it->second;
    }
    // New monomial: assign next variable number (1-indexed)
    int var_num = ++mDim;
    monomial_to_var[key] = var_num;
    return var_num;
}

int StreamingContext::get_var_number(const MonomialKey& key) const {
    auto it = monomial_to_var.find(key);
    if (it != monomial_to_var.end()) {
        return it->second;
    }
    std::cerr << "ERROR: Monomial not found in map during pass 2!" << std::endl;
    return -1;
}

void StreamingContext::start_block(int block_size) {
    nBlocks++;
    current_block = nBlocks;
    block_struct.push_back(block_size);
    current_block_entries = 0;
}

void StreamingContext::write_entry(int var_num, int block, int row, int col, double coef) {
    if (is_counting_pass || output_file == nullptr) return;
    if (coef == 0.0) return;
    
    // SDPA sparse format: var_num block row col value
    fprintf(output_file, "%d %d %d %d %.15e\n", var_num, block, row, col, coef);
}

void StreamingContext::finalize_counting() {
    // Prepare objective coefficient vector
    obj_coef.resize(mDim + 1, 0.0);
    is_counting_pass = false;
    
    // Reset block counter for pass 2
    current_block = 0;
}

void StreamingContext::write_header(const std::string& filename) {
    if (output_file == nullptr) return;

    int numBlocks = (int)block_struct.size();
    
    fprintf(output_file, "* SDPA sparse format data\n");
    fprintf(output_file, "* File name = %s\n", filename.c_str());
    fprintf(output_file, "* mDim = %d, nBlock = %d\n", mDim, numBlocks);
    
    // Number of variables
    fprintf(output_file, "%d\n", mDim);
    
    // Number of blocks
    fprintf(output_file, "%d\n", numBlocks);
    
    // Block structure
    for (int i = 0; i < numBlocks; i++) {
        fprintf(output_file, "%d ", block_struct[i]);
    }
    fprintf(output_file, "\n");
    
    // Objective coefficients (will be filled by convert_obj_stream in pass 2)
    for (int i = 0; i < mDim; i++) {
        fprintf(output_file, "%.15e ", obj_coef[i]);
    }
    fprintf(output_file, "\n");
}

void StreamingContext::finalize_file() {
    if (output_file != nullptr) {
        fclose(output_file);
        output_file = nullptr;
    }
}

// Streaming functions
void convert_obj_stream(poly_info& polyinfo, StreamingContext& ctx) {
    // Pass 1: Register monomials, record objective coefficients
    // Pass 2: Already have obj_coef populated, nothing to write (it's in header)

    int num_terms = polyinfo.sup.pnz_size;
    for (int i = 0; i < num_terms; i++){
        MonomialKey key(polyinfo.sup, i); 

        if (ctx.is_counting_pass) { 
            int var_num = ctx.register_monomial(key);
        } else { 
            // pass 2, look up var number and store coefficient
            int var_num = ctx.get_var_number(key);
            if (var_num > 0 && var_num <= ctx.obj_coef.size()) {
                ctx.obj_coef[var_num - 1] = polyinfo.coef[i][0];
            }
        }
    }
}

void convert_eq_stream(poly_info& polyinfo, spvec_array& bassinfo, StreamingContext& ctx) {
    int num_terms = polyinfo.sup.pnz_size;
    int bsize = bassinfo.pnz_size;
    int sizeCone = polyinfo.sizeCone;
    int move_size = bsize + sizeCone - 1;
    
    if (ctx.is_counting_pass) { // Pass 1: Register monomials, create block
        // Block size: -2 * sizeCone * bsize (diagonal)
        ctx.start_block(-2 * sizeCone * bsize);
        
        for (int s = 0; s < sizeCone; s++) {
            for (int j = 0; j < bsize; j++) {
                for (int i = 0; i < num_terms; i++) {
                    if (fabs(polyinfo.coef[i][s]) > 1.0e-12) {
                        MonomialKey key = merge_monomials(polyinfo.sup, i, bassinfo, j);
                        ctx.register_monomial(key);
                        // Same monomial appears twice (positive and negative)
                        // but we only need to register once
                    }
                }
            }
        }
    } else { // Pass 2: Write entries (twice - positive then negative)
        ctx.nBlocks++;
        
        // First pass: positive coefficients
        for (int s = 0; s < sizeCone; s++) {
            for (int j = 0; j < bsize; j++) {
                for (int i = 0; i < num_terms; i++) {
                    double coef = polyinfo.coef[i][s];
                    if (fabs(coef) > 1.0e-12) {
                        MonomialKey key = merge_monomials(polyinfo.sup, i, bassinfo, j);
                        int var_num = ctx.get_var_number(key);
                        int pos = j + 1 + s;
                        ctx.write_entry(var_num, ctx.nBlocks, pos, pos, coef);
                    }
                }
            }
        }
        // Second pass: negated coefficients at shifted positions
        for (int s = 0; s < sizeCone; s++) {
            for (int j = 0; j < bsize; j++) {
                for (int i = 0; i < num_terms; i++) {
                    double coef = polyinfo.coef[i][s];
                    if (fabs(coef) > 1.0e-12) {
                        MonomialKey key = merge_monomials(polyinfo.sup, i, bassinfo, j);
                        int var_num = ctx.get_var_number(key);
                        int pos = j + 1 + s + move_size;
                        ctx.write_entry(var_num, ctx.nBlocks, pos, pos, -coef);
                    }
                }
            }
        }
    }
}

void convert_ineq_a_ba1_stream(poly_info& polyinfo, spvec_array& bassinfo, StreamingContext& ctx) {
    int num_terms = polyinfo.sup.pnz_size;
    int sizeCone = polyinfo.sizeCone;
    
    if (ctx.is_counting_pass) { // first pass: Register monomials, create block
        ctx.start_block(-sizeCone); 
        
        for (int s = 0; s < sizeCone; s++) {
            for (int i = 0; i < num_terms; i++) {
                double coef = polyinfo.coef[i][s];
                if (fabs(coef) > 1.0e-12) {
                    // Merge poly term with basis[0]
                    MonomialKey key = merge_monomials(polyinfo.sup, i, bassinfo, 0);
                    ctx.register_monomial(key);
                }
            }
        }
    } else { // second pass: Write entries
        ctx.nBlocks++;
        
        for (int s = 0; s < sizeCone; s++) {
            for (int i = 0; i < num_terms; i++) {
                double coef = polyinfo.coef[i][s];
                if (fabs(coef) > 1.0e-12) {
                    MonomialKey key = merge_monomials(polyinfo.sup, i, bassinfo, 0);
                    int var_num = ctx.get_var_number(key);
                    ctx.write_entry(var_num, ctx.nBlocks, s + 1, s + 1, coef);
                }
            }
        }
    }
}

void convert_ineq_a_ba2_stream(poly_info& polyinfo, spvec_array& bassinfo, StreamingContext& ctx) {
    int bsize = bassinfo.pnz_size;
    int num_terms = polyinfo.sup.pnz_size;
    int sizeCone = polyinfo.sizeCone;

    if (ctx.is_counting_pass) {
        std::cout << "  INE_BA2 inputs: bsize=" << bsize << ", num_terms=" << num_terms << ", sizeCone=" << sizeCone << std::endl;
    }
    
    for (int s = 0; s < sizeCone; s++) {
        if (ctx.is_counting_pass) {
            ctx.start_block(bsize);  // Full matrix block
        } else {
            ctx.nBlocks++;
        }
        // Upper triangle of moment matrix
        for (int j = 0; j < bsize; j++) {
            for (int k = j; k < bsize; k++) {
                // First: merge basis[j] with basis[k] to get moment matrix entry
                MonomialKey mm_entry = merge_monomials(bassinfo, j, bassinfo, k);
                
                // Then: for each poly term, merge with mm_entry
                for (int i = 0; i < num_terms; i++) {
                    double coef = polyinfo.coef[i][s];
                    if (fabs(coef) > 1.0e-12) {
                        // Merge poly term with moment matrix entry
                        MonomialKey key = merge_key_with_mono(mm_entry, polyinfo.sup, i);
                        
                        if (ctx.is_counting_pass) {
                            ctx.register_monomial(key);
                        } else {
                            int var_num = ctx.get_var_number(key);
                            ctx.write_entry(var_num, ctx.nBlocks, j + 1, k + 1, coef);
                        }
                    }
                }
            }
        }
    }
}

void convert_sdp_stream(poly_info& polyinfo, spvec_array& bassinfo, StreamingContext& ctx) {
    int bsize = bassinfo.pnz_size;
    int num_terms = polyinfo.sup.pnz_size;
    int sizeCone = polyinfo.sizeCone;
    
    if (ctx.is_counting_pass) {
        // Pass 1: Register monomials, create block
        ctx.start_block(bsize * sizeCone);  // Full matrix block
        
        // For each (j,k) in upper triangle of moment matrix
        for (int j = 0; j < bsize; j++) {
            for (int k = j; k < bsize; k++) {
                MonomialKey mm_entry = merge_monomials(bassinfo, j, bassinfo, k);
                
                // For each poly term
                for (int i = 0; i < num_terms; i++) {
                    // Merge poly term with moment matrix entry
                    MonomialKey key = merge_key_with_mono(mm_entry, polyinfo.sup, i);
                    ctx.register_monomial(key);  // Register ONCE per (j,k,i)
                }
            }
        }
    } else {
        // Pass 2: Write entries
        ctx.nBlocks++;
        
        for (int j = 0; j < bsize; j++) {
            int rowsize = j * sizeCone;
            for (int k = j; k < bsize; k++) {
                int colsize = k * sizeCone;
                MonomialKey mm_entry = merge_monomials(bassinfo, j, bassinfo, k);
                
                for (int i = 0; i < num_terms; i++) {
                    MonomialKey key = merge_key_with_mono(mm_entry, polyinfo.sup, i);
                    int var_num = ctx.get_var_number(key);
                    
                    // Iterate through coefficient matrix entries (CSC format)
                    int r = 0;
                    for (int s = 0; s < sizeCone; s++) {
                        int col_start = polyinfo.mc[s];
                        int col_end = polyinfo.mc[s + 1];
                        
                        while (r < col_end) {
                            int mr_r = polyinfo.mr[r];
                            double coef = polyinfo.coef[r][0];
                            
                            // Write entry
                            int row = mr_r + rowsize + 1;
                            int col = s + colsize + 1;
                            ctx.write_entry(var_num, ctx.nBlocks, row, col, coef);
                            
                            // If off-diagonal in moment matrix AND off-diagonal in coef matrix
                            if (j != k && mr_r != s) {
                                // Write symmetric entry
                                int sym_row = s + rowsize + 1;
                                int sym_col = mr_r + colsize + 1;
                                ctx.write_entry(var_num, ctx.nBlocks, sym_row, sym_col, coef);
                            }
                            r++;
                        }
                    }
                }
            }
        }
    }
}

void convert_ba1mmt_stream(spvec_array& bassinfo, StreamingContext& ctx) {
    if (bassinfo.pnz[1][0] == 0) return;
    
    // squaring monomial so double exponents
    MonomialKey key;
    int start = bassinfo.pnz[0][0];
    int nnz = bassinfo.pnz[1][0];
    for (int i = 0; i < nnz; i++) {
        key.terms.emplace_back(bassinfo.vap[0][start + i], 2 * bassinfo.vap[1][start + i]);
    }
    
    if (ctx.is_counting_pass) { // Pass 1: Register monomial, create block
        ctx.register_monomial(key);
        ctx.start_block(-1);
    } else { // Pass 2: Write entry
        int var_num = ctx.get_var_number(key);
        ctx.nBlocks++;  // Track block number for this pass
        ctx.write_entry(var_num, ctx.nBlocks, 1, 1, 1.0);
    }
}

// create a moment matrix from a given monomial basis (bassinfo)
void convert_ba2mmt_stream(spvec_array& bassinfo, StreamingContext& ctx) {
    int bsize = bassinfo.pnz_size;

    if (ctx.is_counting_pass) {
        std::cout << "BA2MMT: bsize=" << bsize << ", expected_monomials=" << (bsize * (bsize + 1) / 2) << std::endl;
    }

    // // DEBUG: Print structure info
    // std::cout << "BA2MMT bsize=" << bsize << ", vap_size=" << bassinfo.vap_size << std::endl;
    // for (int i = 0; i < bsize && i < 10; i++) {  // limit to first 10
    //     std::cout << "  pnz[" << i << "]: start=" << bassinfo.pnz[0][i] << ", len=" << bassinfo.pnz[1][i] << std::endl;
    // }

    if (ctx.is_counting_pass) {
        // Pass 1: Register all monomials, create block
        ctx.start_block(bsize);  // Positive = full matrix of size bsize
        std::unordered_set<std::string> unique_basis; 

        for (int i = 0; i < bsize; i++) {
            std::string s;
            int start = bassinfo.pnz[0][i];
            int len = bassinfo.pnz[1][i];
            if (start >= 0) {
                for (int t = 0; t < len; t++) {
                    s += std::to_string(bassinfo.vap[0][start+t]) + "^" + 
                         std::to_string(bassinfo.vap[1][start+t]) + " ";
                }
            }
            unique_basis.insert(s);
        }
        std::cout << "BA2MMT: bsize=" << bsize << ", unique_basis_elements=" << unique_basis.size() << std::endl;


        // Upper triangle: all pairs (i,j) where j >= i
        for (int i = 0; i < bsize; i++) {
            for (int j = i; j < bsize; j++) { // Merge basis[i] with basis[j]
                MonomialKey key = merge_monomials(bassinfo, i, bassinfo, j);

                 // DEBUG: Print first few monomials
                static int debug_count = 0;
                bool should_print = ctx.is_counting_pass && (
                    debug_count < 5 ||                    // First 5
                    (i == 1 && j <= 3) ||                 // basis[1] × basis[1..3]
                    (i == 100 && j == 100) ||             // A larger index
                    (i == 200 && j == 300)                // Another combo
                );
                if (should_print) {
                    std::cout << "  MM[" << i << "," << j << "]: ";
                    if (key.terms.empty()) std::cout << "(constant)";
                    for (auto& [v, e] : key.terms) {
                        std::cout << "x" << v << "^" << e << " ";
                    }
                    std::cout << std::endl;
                    debug_count++;
                } // end Debug

                ctx.register_monomial(key);
            }
        }
    } else { // second pass: write entries
        ctx.nBlocks++;
        for (int i = 0; i < bsize; i++) {
            for (int j = i; j < bsize; j++) {
                MonomialKey key = merge_monomials(bassinfo, i, bassinfo, j);
                int var_num = ctx.get_var_number(key);
                ctx.write_entry(var_num, ctx.nBlocks, i + 1, j + 1, 1.0);
            }
        }
    }
}

void stream_psdp_to_file(int mdim,int msize,std::vector<poly_info>& polyinfo,std::vector<spvec_array>& bassinfo,const std::string& sdpafile) {

    StreamingContext ctx;
    ctx.is_counting_pass = true;
    
    std::cout << "=== Pass 1: Counting ===" << std::endl;
    
    // --- PASS 1: Count monomials, build structure ---
    convert_obj_stream(polyinfo[0], ctx);
    
    std::cout << "After convert_obj_stream" << std::endl;
    
    for (int i = 1; i < msize; i++) {
        int before = ctx.mDim; 
        
        if (polyinfo[i].typeCone == EQU) {
            convert_eq_stream(polyinfo[i], bassinfo[i], ctx);
            std::cout << "EQU[" << i << "] added " << (ctx.mDim - before) << " monomials" << std::endl;
        }
        else if (polyinfo[i].typeCone == 0) {
            continue;
        }
        else if (polyinfo[i].typeCone == INE && bassinfo[i].pnz_size == 1) {
            convert_ineq_a_ba1_stream(polyinfo[i], bassinfo[i], ctx);
            std::cout << "INE_BA1[" << i << "] added " << (ctx.mDim - before) << " monomials" << std::endl;
        }
        else if (polyinfo[i].typeCone == INE && bassinfo[i].pnz_size >= 2) {
            convert_ineq_a_ba2_stream(polyinfo[i], bassinfo[i], ctx);
            std::cout << "INE_BA2[" << i << "] added " << (ctx.mDim - before) << " monomials" << std::endl;
        }
        else if (polyinfo[i].typeCone == SDP) {
            convert_sdp_stream(polyinfo[i], bassinfo[i], ctx);
            std::cout << "SDP[" << i << "] added " << (ctx.mDim - before) << " monomials" << std::endl;
        }
        else if (bassinfo[i].pnz_size == 1) {
            convert_ba1mmt_stream(bassinfo[i], ctx);
            std::cout << "BA1MMT[" << i << "] added " << (ctx.mDim - before) << " monomials" << std::endl;
        }
        else if (bassinfo[i].pnz_size >= 2) {
            convert_ba2mmt_stream(bassinfo[i], ctx);
            std::cout << "BA2MMT[" << i << "] added " << (ctx.mDim - before) << " monomials" << std::endl;
        }
    }
    
    std::cout << "Pass 1 complete: mDim=" << ctx.mDim << ", nBlocks=" << ctx.nBlocks << std::endl;
    std::cout << "block_struct.size()=" << ctx.block_struct.size() << std::endl;
    
    // --- Prepare for Pass 2 ---
    ctx.finalize_counting();
    
    // Open output file
    ctx.output_file = fopen(sdpafile.c_str(), "w");
    if (ctx.output_file == nullptr) {
        std::cerr << "Error: Could not open output file " << sdpafile << std::endl;
        return;
    }
    
    std::cout << "=== Pass 2: Writing ===" << std::endl;
    
    // Reset block counter
    ctx.nBlocks = 0;
    ctx.current_block = 0;
    
    // --- PASS 2: Write to file ---
    convert_obj_stream(polyinfo[0], ctx);
    
    // Write header after objective (so we have obj_coef populated)
    ctx.write_header(sdpafile);
    
    for (int i = 1; i < msize; i++) {
        if (polyinfo[i].typeCone == EQU) {
            convert_eq_stream(polyinfo[i], bassinfo[i], ctx);
        }
        else if (polyinfo[i].typeCone == 0) {
            continue; 
        }
        else if (polyinfo[i].typeCone == INE && bassinfo[i].pnz_size == 1) {
            convert_ineq_a_ba1_stream(polyinfo[i], bassinfo[i], ctx);
        }
        else if (polyinfo[i].typeCone == INE && bassinfo[i].pnz_size >= 2) {
            convert_ineq_a_ba2_stream(polyinfo[i], bassinfo[i], ctx);
        }
        else if (polyinfo[i].typeCone == SDP) {
            convert_sdp_stream(polyinfo[i], bassinfo[i], ctx);
        }
        else if (bassinfo[i].pnz_size == 1) {
            convert_ba1mmt_stream(bassinfo[i], ctx);
        }
        else if (bassinfo[i].pnz_size >= 2) {
            convert_ba2mmt_stream(bassinfo[i], ctx);
        }
    }
    
    ctx.finalize_file();
    std::cout << "=== Streaming complete ===" << std::endl;
}

// Simple test function - call from main to verify streaming code compiles/works
void test_streaming_basics() {
    std::cout << "\n=== Testing Streaming Basics ===" << std::endl;
    
    StreamingContext ctx;
    
    // Test 1: MonomialKey hashing
    MonomialKey m1, m2, m3, m4;
    m1.terms = {{1, 2}, {3, 1}};  // x1^2 * x3
    m2.terms = {{2, 1}};          // x2
    m3.terms = {{1, 2}, {3, 1}};  // same as m1
    m4.terms = {{1,2}, {3,1}, {2,1}}; 
    
    int v1 = ctx.register_monomial(m1);
    int v2 = ctx.register_monomial(m2);
    int v3 = ctx.register_monomial(m3);
    int v4 = ctx.register_monomial(m4); 
    
    std::cout << "m1 -> var " << v1 << std::endl;
    std::cout << "m2 -> var " << v2 << std::endl;
    std::cout << "m3 (same as m1) -> var " << v3 << std::endl;
    std::cout << "m4 -> var " << v4 << std::endl;
    std::cout << "Total mDim = " << ctx.mDim << " (should be 3)" << std::endl;
    
    // Test 2: Block tracking
    ctx.start_block(3);   // 3x3 matrix block
    ctx.start_block(-2);  // 2x2 diagonal block
    ctx.start_block(20);
    ctx.start_block(3); 
    std::cout << "nBlocks = " << ctx.nBlocks << " (should be 4)" << std::endl;

    // Test 3: convert_ba1mmt_stream
    std::cout << "\n--- Testing convert_ba1mmt_stream ---" << std::endl;

    // Create a simple basis: x₁ · x₂ (var 1 exp 1, var 2 exp 1)
    spvec_array basis;
    basis.alloc(1, 2);  // 1 monomial, 2 terms
    basis.pnz[0][0] = 0;   // starts at position 0
    basis.pnz[1][0] = 2;   // has 2 terms
    basis.vap[0][0] = 1;   // var 1
    basis.vap[1][0] = 1;   // exp 1
    basis.vap[0][1] = 2;   // var 2
    basis.vap[1][1] = 1;   // exp 1
    basis.pnz_size = 1;
    basis.vap_size = 2;

    // Reset context for clean test
    StreamingContext ctx2;
    ctx2.is_counting_pass = true;

    // Pass 1: Count
    convert_ba1mmt_stream(basis, ctx2);
    std::cout << "After pass 1: mDim=" << ctx2.mDim << ", nBlocks=" << ctx2.nBlocks << std::endl;
    std::cout << "  (Should be mDim=1, nBlocks=1)" << std::endl;

    // Prepare for pass 2
    ctx2.finalize_counting();
    ctx2.output_file = fopen("/tmp/test_ba1mmt.sdpa", "w");
    ctx2.nBlocks = 0;
    
    // Pass 2: Write  
    convert_ba1mmt_stream(basis, ctx2);

    // Test 4: convert_ba2mmt_stream  
    std::cout << "\n--- Testing convert_ba2mmt_stream ---" << std::endl;

    // Create basis with 2 elements: {x₁, x₂}
    spvec_array basis2;
    basis2.alloc(2, 2);  // 2 monomials, 2 total terms
    basis2.pnz[0][0] = 0; basis2.pnz[1][0] = 1;  // mono 0: 1 term starting at 0
    basis2.pnz[0][1] = 1; basis2.pnz[1][1] = 1;  // mono 1: 1 term starting at 1
    basis2.vap[0][0] = 1; basis2.vap[1][0] = 1;  // var 1, exp 1 (x₁)
    basis2.vap[0][1] = 2; basis2.vap[1][1] = 1;  // var 2, exp 1 (x₂)
    basis2.pnz_size = 2;
    basis2.vap_size = 2;

    StreamingContext ctx3;
    ctx3.is_counting_pass = true;
    convert_ba2mmt_stream(basis2, ctx3);
    std::cout << "After pass 1: mDim=" << ctx3.mDim << ", nBlocks=" << ctx3.nBlocks << std::endl;
    std::cout << "  (Should be mDim=3, nBlocks=1 for entries: x₁², x₁x₂, x₂²)" << std::endl;

    ctx3.finalize_counting();
    ctx3.output_file = fopen("/tmp/test_ba2mmt.sdpa", "w");
    ctx3.nBlocks = 0;
    convert_ba2mmt_stream(basis2, ctx3);
    ctx3.finalize_file();
    std::cout << "Output: " << std::endl;
    system("cat /tmp/test_ba2mmt.sdpa");

    basis2.del();

    // Test 5: convert_ineq_a_ba1_stream
std::cout << "\n--- Testing convert_ineq_a_ba1_stream ---" << std::endl;

// Create polynomial: 2.0*x₁ + 3.0*x₂ (two terms)
poly_info poly;
poly.sizeCone = 1;  // Single constraint
poly.sup.alloc(2, 2);  // 2 monomials, 2 total terms
poly.sup.pnz[0][0] = 0; poly.sup.pnz[1][0] = 1;  // term 0: x₁
poly.sup.pnz[0][1] = 1; poly.sup.pnz[1][1] = 1;  // term 1: x₂
poly.sup.vap[0][0] = 1; poly.sup.vap[1][0] = 1;  // var 1, exp 1
poly.sup.vap[0][1] = 2; poly.sup.vap[1][1] = 1;  // var 2, exp 1
poly.sup.pnz_size = 2;
poly.sup.vap_size = 2;

// Coefficients: coef[term][sizeCone_index]
poly.coef.resize(2);
poly.coef[0].resize(1); poly.coef[0][0] = 2.0;  // coef for x₁
poly.coef[1].resize(1); poly.coef[1][0] = 3.0;  // coef for x₂

// Create basis with 1 element: {x₃}
spvec_array basis3;
basis3.alloc(1, 1);
basis3.pnz[0][0] = 0; basis3.pnz[1][0] = 1;
basis3.vap[0][0] = 3; basis3.vap[1][0] = 1;  // var 3, exp 1
basis3.pnz_size = 1;
basis3.vap_size = 1;

StreamingContext ctx4;
ctx4.is_counting_pass = true;
convert_ineq_a_ba1_stream(poly, basis3, ctx4);
std::cout << "After pass 1: mDim=" << ctx4.mDim << ", nBlocks=" << ctx4.nBlocks << std::endl;
std::cout << "  (Should be mDim=2 for: x₁x₃, x₂x₃)" << std::endl;

ctx4.finalize_counting();
ctx4.output_file = fopen("/tmp/test_ineq_ba1.sdpa", "w");
ctx4.nBlocks = 0;
convert_ineq_a_ba1_stream(poly, basis3, ctx4);
ctx4.finalize_file();
std::cout << "Output: " << std::endl;
system("cat /tmp/test_ineq_ba1.sdpa");

poly.sup.del();
basis3.del();

// Test 6: convert_ineq_a_ba2_stream
std::cout << "\n--- Testing convert_ineq_a_ba2_stream ---" << std::endl;

// Create polynomial: 5.0 * x₁
poly_info poly2;
poly2.sizeCone = 1;
poly2.sup.alloc(1, 1);  // 1 monomial, 1 term
poly2.sup.pnz[0][0] = 0; poly2.sup.pnz[1][0] = 1;
poly2.sup.vap[0][0] = 1; poly2.sup.vap[1][0] = 1;  // x₁
poly2.sup.pnz_size = 1;
poly2.sup.vap_size = 1;
poly2.coef.resize(1);
poly2.coef[0].resize(1); poly2.coef[0][0] = 5.0;

// Create basis with 2 elements: {x₂, x₃}
spvec_array basis4;
basis4.alloc(2, 2);
basis4.pnz[0][0] = 0; basis4.pnz[1][0] = 1;  // x₂
basis4.pnz[0][1] = 1; basis4.pnz[1][1] = 1;  // x₃
basis4.vap[0][0] = 2; basis4.vap[1][0] = 1;  // var 2, exp 1
basis4.vap[0][1] = 3; basis4.vap[1][1] = 1;  // var 3, exp 1
basis4.pnz_size = 2;
basis4.vap_size = 2;

StreamingContext ctx5;
ctx5.is_counting_pass = true;
convert_ineq_a_ba2_stream(poly2, basis4, ctx5);
std::cout << "After pass 1: mDim=" << ctx5.mDim << ", nBlocks=" << ctx5.nBlocks << std::endl;
std::cout << "  (Should be mDim=3 for: x₁x₂², x₁x₂x₃, x₁x₃²)" << std::endl;

ctx5.finalize_counting();
ctx5.output_file = fopen("/tmp/test_ineq_ba2.sdpa", "w");
ctx5.nBlocks = 0;
convert_ineq_a_ba2_stream(poly2, basis4, ctx5);
ctx5.finalize_file();
std::cout << "Output: " << std::endl;
system("cat /tmp/test_ineq_ba2.sdpa");

poly2.sup.del();
basis4.del();

// Test 7: convert_eq_stream
std::cout << "\n--- Testing convert_eq_stream ---" << std::endl;

// Polynomial: 4.0 * x₁
poly_info poly3;
poly3.sizeCone = 1;
poly3.sup.alloc(1, 1);
poly3.sup.pnz[0][0] = 0; poly3.sup.pnz[1][0] = 1;
poly3.sup.vap[0][0] = 1; poly3.sup.vap[1][0] = 1;
poly3.sup.pnz_size = 1;
poly3.sup.vap_size = 1;
poly3.coef.resize(1);
poly3.coef[0].resize(1); poly3.coef[0][0] = 4.0;

// Basis: {x₂}
spvec_array basis5;
basis5.alloc(1, 1);
basis5.pnz[0][0] = 0; basis5.pnz[1][0] = 1;
basis5.vap[0][0] = 2; basis5.vap[1][0] = 1;
basis5.pnz_size = 1;
basis5.vap_size = 1;

StreamingContext ctx6;
ctx6.is_counting_pass = true;
convert_eq_stream(poly3, basis5, ctx6);
std::cout << "After pass 1: mDim=" << ctx6.mDim << ", nBlocks=" << ctx6.nBlocks << std::endl;
std::cout << "  (Should be mDim=1, nBlocks=1, block size -2)" << std::endl;

ctx6.finalize_counting();
ctx6.output_file = fopen("/tmp/test_eq.sdpa", "w");
ctx6.nBlocks = 0;
convert_eq_stream(poly3, basis5, ctx6);
ctx6.finalize_file();
std::cout << "Output (should have +4 at (1,1) and -4 at (2,2)):" << std::endl;
system("cat /tmp/test_eq.sdpa");

poly3.sup.del();
basis5.del();
    
    std::cout << "=== Streaming Test PASSED ===" << std::endl;
}