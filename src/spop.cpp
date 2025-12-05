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
    
    // Solve with SDPA
    if(POP.param.SDPsolverSW == 1 && info.infeasibleSW == 0) {
        SDPA Problem;
        MakeSDPAform(sdpdata, Problem);
        
        Problem.setDisplay(NULL);
        Problem.initializeSolve();
        Problem.setParameterEpsilonStar(POP.param.SDPsolverEpsilon);
        Problem.setParameterEpsilonDash(POP.param.SDPsolverEpsilon);
        
        Problem.solve();
        Problem.terminate();
    }
}