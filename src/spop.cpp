#include "spop.h"

#include "global.h"
#include "conversion.h"
#include "polynomials.h"
#include "sup.h"
#include "spvec.h"
#include "input.h"
#include "info.h"
#include "cspGraph.h"
#include "debug.h"
#include "Parameters.h"
#include "sdpa_call.h"

#include <iostream>
#include <fstream>

extern void makeSDPr(s3r& POP, mysdp& sdpdata, Info& info, std::string& pname, std::vector<std::vector<double>>& fixedVar, std::tuple<int,int, std::vector<int>, std::vector<int>, std::vector<int>>& fromGen);
extern void MakeSDPAform(mysdp& sdpdata, SDPA& Problem);

void solveWithSparsePOP(std::string& gmsFilePath, std::tuple<int,int, std::vector<int>, std::vector<int>, std::vector<int>>& fromGen) {
    // Create SparsePOP objects
    s3r POP;
    mysdp sdpdata;
    Info info;
    
    // Convert POP to SDP
    std::vector<std::vector<double>> fixedVar(2);
    makeSDPr(POP, sdpdata, info, gmsFilePath, fixedVar, fromGen);

    std::cout << "\n=== SDP Problem Info ===" << std::endl;
    std::cout << "SDP dimensions: mDim = " << sdpdata.mDim << ", nBlocks = " << sdpdata.nBlocks << std::endl;
    std::cout << "- Number of variables: " << POP.Polysys.dimVar << std::endl;
    std::cout << "- Number of constraints: " << POP.Polysys.numSys << std::endl;
    std::cout << "nBlocks = " << sdpdata.nBlocks << std::endl;
    for(int i=1; i <= sdpdata.nBlocks; i++){
        if(sdpdata.block_info[1][i] > 0){
            std::cout << "Block " << i << ": size = " << abs(sdpdata.bLOCKsTruct[i]) << std::endl;
        }
    }
    if (info.infeasibleSW != 0) {
        std::cout << "WARNING: Problem detected as infeasible before solving!" << std::endl;
        return;
    }
    std::cout << "========================\n" << std::endl;
    
    // Solve with SDPA
    if(POP.param.SDPsolverSW == 1 && info.infeasibleSW == 0) {
        SDPA Problem;
        std::cout << "Before calling MakeSDPAform " << std::endl;
        MakeSDPAform(sdpdata, Problem);
        std::cout << "After calling MakeSDPAform " << std::endl;

        Problem.setDisplay(NULL);
        Problem.initializeSolve();
        Problem.setParameterEpsilonStar(POP.param.SDPsolverEpsilon);
        Problem.setParameterEpsilonDash(POP.param.SDPsolverEpsilon);
        
        Problem.solve();
        Problem.terminate();
    }
}