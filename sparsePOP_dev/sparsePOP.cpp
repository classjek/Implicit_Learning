/* -------------------------------------------------------------
 *
 * This file is a component of SparsePOP
 * Copyright (C) 2007 SparsePOP Project
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 *
 * ------------------------------------------------------------- */
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
//#include "checkLinDep.h"

#include "sdpa_call.h"

/* convert our form into SDPA sparse format */
void MakeSDPAform(class mysdp & sdpdata, SDPA & Problem){
    
    //cout << "test" << endl;
    char pname[11];
    strcpy(pname, "param.sdpa");
    char *pname_ptr = &(pname[0]);
    //printf("%s\n", pname);
    //FILE* fpresult;
    //fpresult=fopen("test.out", "r");
    Problem.readParameter(pname_ptr);
    //strcpy(Problem.ParameterFileName,"param.sdpa");
    //Problem.ParameterFile = fopen(Problem.ParameterFileName,"r");
    //cout << "1" <<endl;
    //cout << "2" <<endl;
    //cout << "3" <<endl;
    Problem.inputConstraintNumber(sdpdata.mDim);
    Problem.inputBlockNumber(sdpdata.nBlocks);
    //cout << "mDim   = " << sdpdata.mDim << endl;
    //cout << "nBlock = " << Problem.nBLOCK << endl;
    for(int i=0; i<sdpdata.nBlocks; i++){
        //cout << "i= ";
        //cout << sdpdata.bLOCKsTruct[i+1] << endl;
        if(sdpdata.bLOCKsTruct[i+1] > 0){
            Problem.inputBlockSize(i+1, sdpdata.bLOCKsTruct[i+1]);
            Problem.inputBlockType(i+1, SDPA::SDP);
        }else{
            Problem.inputBlockSize(i+1, -sdpdata.bLOCKsTruct[i+1]);
            Problem.inputBlockType(i+1, SDPA::LP);
        }
    }
    //cout << "4" <<endl;
    //cout << "5" <<endl;
    Problem.initializeUpperTriangleSpace();
    for(int j=sdpdata.block_info[0][0];j<sdpdata.block_info[0][0]+sdpdata.block_info[1][0];j++){
        Problem.inputCVec(sdpdata.ele.sup.pnz[0][j], sdpdata.ele.coef[j]);
        //SDPA_Input_cVECT(Problem,sdpdata.ele.sup.pnz[0][j],sdpdata.ele.coef[j]);
    }
    //cout << "6" <<endl;
    vector<vector<int> > nnz(1+sdpdata.mDim);
    for(int i=0; i<1+sdpdata.mDim; i++){
        nnz[i].resize(sdpdata.nBlocks);
        for(int j=0; j<sdpdata.nBlocks; j++){
            nnz[i][j] = 0;
        }
    }
    //cout << "7" <<endl;
    /*
     * cout << "pnz = " << endl;
     * for(int i=0; i<sdpdata.ele.sup.pnz[0].size(); i++){
     * cout << "pnz[0][" << i << "]= " << sdpdata.ele.sup.pnz[0][i] << endl;
     * cout << "pnz[1][" << i << "]= " << sdpdata.ele.sup.pnz[1][i] << endl;
     * }
     */
    int row1, col1;
    for(int i=1;i<sdpdata.nBlocks+1;i++){
        for(int j=sdpdata.block_info[0][i];j<sdpdata.block_info[0][i]+sdpdata.block_info[1][i];j++){
            row1 = sdpdata.ele.sup.pnz[0][j];
            col1 = sdpdata.ele.bij[0][j]-1;
            nnz[row1][col1]++;
            //cout << "nnz[" << row1 << "][" << col1 << "]= " << nnz[row1][col1] << endl;
        }
    }
    //cout << "8" <<endl;
    /*
     * cout << "nnz information " << endl;
     * for(int i=0; i<1+Problem.mDIM; i++){
     * for(int j=0;j<Problem.nBLOCK;j++){
     * cout << "nnz[" << i << "][" << j << "] = " << nnz[i][j] << " ";
     * }
     * cout << endl;
     * }
     */
    for(int i=0; i<1+sdpdata.mDim;i++){
        for(int j=0;j<sdpdata.nBlocks;j++){
            //SDPA_CountUpperTriangle(Problem,i,j+1,nnz[i][j]);
        }
    }
    //cout << "9" <<endl;
    //system("top -b -n 1|grep sparsePOP |head -1|awk '{printf(\"memory = %s\\n\"),$6}' ");
    //SDPA_Make_sfMAT(Problem);
    //cout << "10" <<endl;
    int no;
    int bn;
    int row;
    int col;
    double coef;
    for(int i=1;i<sdpdata.nBlocks+1;i++){
        for(int j=sdpdata.block_info[0][i];j<sdpdata.block_info[0][i]+sdpdata.block_info[1][i];j++){
            no  = sdpdata.ele.sup.pnz[0][j];
            bn  = sdpdata.ele.bij[0][j];
            row = sdpdata.ele.bij[1][j];
            col = sdpdata.ele.bij[2][j];
            coef= sdpdata.ele.coef[j];
            //SDPA_InputElement(Problem,no,bn,row,col,coef);
            Problem.inputElement(no, bn, row, col, coef);
        }
    }
    Problem.initializeUpperTriangle();
    Problem.setInitPoint(false);
	Problem.setNumThreads(0);	
	cout << "## SDPA uses " << Problem.getNumThreads() << " threads." <<endl;
    //cout << "11" <<endl;
    /*
     * char oname[9];
     * strcpy(oname, "test.out");
     * char *oname_ptr = &(oname[0]);
     * printf("%s\n", pname);
     * char oform[6];
     * strcpy(oform,"%5.3e");
     * char *oform_ptr = &(oform[0]);
     * Problem.writeInputSparse(oname_ptr, oform_ptr);
     */
}

/* Input results of SDPA */
void InputResultSDPA(SDPA & Problem, class Info & info, double ObjscaleVal, double ObjConst){
    /*
     * cout.precision(15);
     * cout << "PrimalObj   = " << Problem.PrimalObj << endl;
     * cout << "ObjscaleVal = " << ObjscaleVal << endl;
     * cout << "ObjConst    = " << ObjConst << endl;
     */
    info.PrimalObj   = Problem.getPrimalObj();
    info.DualObj     = Problem.getDualObj();
    info.SDPobjValue = info.PrimalObj*ObjscaleVal + ObjConst;
    //info.SDPobjValue = Problem.PrimalObj*ObjscaleVal + ObjConst;
    info.PrimalError = Problem.getPrimalError();
    info.DualError   = Problem.getDualError();
    info.SDPAtime    = Problem.getSolveTime();
    info.iter        = Problem.getIteration();
    info.digit       = Problem.getDigits();
    char phase_string[30];
    Problem.getPhaseString(phase_string);
    info.sdpaPhase   = phase_string;
    double m = fabs(info.PrimalObj) + fabs(info.DualObj);
    m = m/2;
    if(m > 1){
        m = 1;
    }
    info.relativeGap = fabs(info.PrimalObj-info.DualObj)/m;
    info.mu = Problem.getMu();
    int n =0;
    for(int i=0; i<Problem.getBlockNumber(); i++){
        if(Problem.getBlockSize(i+1) > 0){
            n = n + Problem.getBlockSize(i+1);
        }else{
            n = n - Problem.getBlockSize(i+1);
        }
    }
    info.gap = info.mu*n;
}

void makeSDPr(class s3r &POP, class mysdp & sdpdata, class Info & info, string & pname, vector<vector<double> > & fixedVar){
    /* read GMS file */
    POP.problemName = pname;

    // testing delete me later
    tuple<int,int, vector<int>, vector<int>, vector<int>> fromGen = loadBlob();

    inputGMS(POP.Polysys, pname);
    //cout << "reading input finished. " << endl;
    /* read param file */
    POP.param.SetParameters("param.pop", POP.Polysys.dimVar);

    int maxD = POP.Polysys.polynomial[0].degree;
    for(int i=1; i<POP.Polysys.numSys; i++){
        if(maxD < POP.Polysys.polynomial[i].degree){
            maxD = POP.Polysys.polynomial[i].degree;
        }
    }
    double maxD2 = maxD/2;
    if(POP.param.relax_Order <maxD2){
        if(maxD % 2 == 0){
            POP.param.relax_Order = (maxD >> 1);
        }else{
            POP.param.relax_Order = ((maxD+1) >> 1);
        }
    }
    
    /* allocate for POPs */
	if(POP.Polysys.minOrmax == -1){
		cout << "# The inputted problem is the maximization problem." << endl;
		cout << "# We multiply -1 into the objective fuction for" << endl; 
		cout << "# converting it into the minimization problem." << endl;
	}
    POP.objPoly = POP.Polysys.polynomial[0];
    POP.detailedInfFile = POP.param.detailedInfFile;
    POP.sdpaDataFile = POP.param.sdpaDataFile;
    
    /* save original POP */
    POP.OriPolysys = POP.Polysys;
    //POP.OriPolysys.writePolynomials();
    
    int ndim = POP.Polysys.dimVar;
    fixedVar[0].resize(ndim, 0);
    fixedVar[1].resize(ndim, 10*MIN);

    vector<set<int>> expectMap(POP.Polysys.dimVar);
    //deleteVar(POP.Polysys, fixedVar, varMapping);
	
    POP.degOneTerms.resize(POP.Polysys.dimVar, 0);
    if(POP.Polysys.dimVar == 0){
        //cout << POP.Polysys.dimVar << endl;
        info.nDimZero(POP.OriPolysys, fixedVar);
        vector<double> fValue;
        vector<double> maxAbs;
        POP.OriPolysys.polynomial[0].evalPoly(info.xvect, fValue, maxAbs);
        info.objValue = fValue[0];
        if(info.infeasibleSW == 2 || info.infeasibleSW == -2){
            return;
        }
    }
    double objconst;
    int slen = POP.Polysys.numSys;
    int blen = POP.Polysys.dimVar;
    int mlen = POP.Polysys.dimVar * POP.Polysys.dimVar;
    vector<double> scalevalue(slen);
    vector<double> bvect(blen);
    vector<double> permmatrix(mlen);
    
    conversion_part1(POP, objconst,
            slen, scalevalue,
            blen, bvect,
            mlen, permmatrix);
    //cout << "conversion1 finished. " << endl;
    //POP.Polysys.writePolynomials();
    POP.scalevalue = scalevalue;
    POP.bvect = bvect;
    POP.permmatrix =permmatrix;
    //printScaleInfo(slen,blen,mlen,scalevalue,bvect,permmatrix);
    vector<vector<double> > cspmat;
    genCSP(POP.Polysys, cspmat);
	//printCsp(cspmat);
    
    vector<int> oriidx(POP.Polysys.dimVar);
    ///vector<vector<int> > extmat(POP.Polysys.dimVar);
    class SparseMat extmat;
    
    if(POP.param.sparseSW == 1){
        if(POP.param.Method == "metis"){
            MetisGraph graph;
            graph.makeGraph(cspmat);
            METIS_NodeND(&graph.ndim, graph.xadj, graph.adjncy, &graph.numflag, graph.options, graph.perm, graph.iperm);
            for(int i=0; i< POP.Polysys.dimVar; i++){
                oriidx[i] = graph.perm[i];
            }
        }else if(POP.param.Method == "symamd"){
            SymamdGraph graph;
            graph.makeGraph(cspmat);
            int ok = symamd(graph.ndim, graph.B, graph.q, graph.perm, (double *)NULL, graph.stats, &calloc, &free);
            //symamd_report(graph.stats);
			if(ok != 1){
				printf("symamd fliled. Construct a postive definite csp matrix.");
				exit(EXIT_FAILURE);
			}	
            for(int i=0; i<graph.ndim; i++){
                oriidx[i] = graph.perm[i];
		//cout << "oriidx[" << i <<"]=" << oriidx[i] <<endl;
            }
            //graph.disp();
        }else{
            cout << " ## Error: Should be metis or symamd in param.Method." << endl;
            exit(EXIT_FAILURE);
        }
        char UPLO = 'L';
        int info;
        int size = POP.Polysys.dimVar*(POP.Polysys.dimVar+1) >> 1;
        double *A = new double[size];
        makeA(A, cspmat, oriidx);
        //printMat(A,POP.Polysys.dimVar);
        
        /* Cholesky decomposition */
        /* Fiding chordal extension */
        dpptrf_(&UPLO, &POP.Polysys.dimVar, A, &info);
        if(info == 0){
            //printMat(A,POP.Polysys.dimVar);
            makeExtmat(POP.Polysys.dimVar, A, extmat);
        }else if(info <0){
            cout << " ## Error: csp matrix has something wrong." << endl;
            exit(EXIT_FAILURE);
        }else if(info > 0){
            cout << " ## Error: csp matrix is not positive definite." << endl;
            cout << "           Try to cosntruct csp matrix again." << endl;
            exit(EXIT_FAILURE);
        }
        delete [] A;
        A = NULL;
    }else if(POP.param.sparseSW == 0){
	int nnz = POP.Polysys.dimVar^2;
	extmat.resizeIr(nnz, 0);
	for(int i=0; i<nnz; i++){
		extmat.ir[i] = i % POP.Polysys.dimVar;
	}
	extmat.resizeJc(POP.Polysys.dimVar+1, 0);
        for(int i=0; i< POP.Polysys.dimVar; i++){
            oriidx[i] = i;
		extmat.jc[i] = i*(POP.Polysys.dimVar);
        }
	extmat.jc[POP.Polysys.dimVar] = nnz;
    }else{
        cout << " ## Error: sparseSW should be 1 or 0. " << endl;
        exit(EXIT_FAILURE);
    }
    
    /* allocate sdpdta */
    // add fixedVar here so we cna acces information regarding which variables
    // were transformed into which other variables during optimization process
    conversion_part2(POP, fixedVar, expectMap, fromGen, oriidx, extmat, sdpdata);
}

// Comment out main for compilation with our code 
// int main(int argc, char *argv[]){
    
//     printf("\nSparsePOP by H.Waki, S.Kim, M.Kojima, M.Muramatsu\n");
//     printf("                       H.Sugimoto and M.Yamashita\n");
//     printf("                                    February 2012\n\n");
//     //printf("                                      August 2011\n\n");
//     //printf("                                        July 2011\n\n");
//     //printf("                                        June 2011\n\n");
//     //printf("                                     March 2011\n\n");
//     //printf("                                          October 2010\n\n");
    
//     if(argc != 2){
//         cout << " ## Error: should input one file name. " << endl;
//         exit(EXIT_FAILURE);
//     }
//     string pname = argv[1];
//     /* file check */
//     ifstream fin;
//     fin.open(argv[1]);
//     if(!fin){
//         cout << " ## Error: cannno open file(" << pname << ")" << endl;
//         exit(EXIT_FAILURE);
//     }else{
//         //cout << pname << endl;
//     }
//     fin.close();
//     class s3r POP;
//     class mysdp sdpdata;
//     class Info info;

    
//     double cputime1, cputime2, cputime3;
//     cputime1 = (double)clock()/(double)CLOCKS_PER_SEC;
//     vector<vector<double> > fixedVar(2);
// 	#ifdef DEBUG
// 	    cout << "making SDPr started " << endl;
// 	#endif
//     makeSDPr(POP, sdpdata, info, pname, fixedVar);
//     if(POP.param.sdpaDataFile.empty() == false){
//         write_sdpa(/*IN*/sdpdata,/*OUT*/ POP.param.sdpaDataFile);
//     }
// 	#ifdef DEBUG
//     cout << "making SDPr finished " << endl;
// 	#endif
//     cputime2 = (double)clock()/(double)CLOCKS_PER_SEC;
//     info.conversion = cputime2-cputime1;
//     //cout << "reduceMomentMatSW_= " << POP.param.reduceMomentMatSW << endl; // testing

//     //sdpdata.disp();
//     if(info.infeasibleSW == 2||info.infeasibleSW == -2){
//         cputime3 = (double)clock()/(double)CLOCKS_PER_SEC;
//         info.total = cputime3 - cputime1;
//         vector<double> fValue;
//         vector<double> maxAbs;
//         POP.OriPolysys.polynomial[0].evalPoly(info.xvect, fValue, maxAbs);
//         info.objValue = fValue[0];
//         info.printSolutions(POP, POP.param.printFileName);
//         if(POP.param.printFileName.empty() == false && POP.param.printOnScreen == 1){
//             info.printSolutions(POP, "");
//         }
//         return 0;
//     }
//     FILE* fpresult = NULL;
//     if(POP.param.SDPsolverSW == 1 && info.infeasibleSW == 0){
//         /* Solve SDP by SDPA */
//         SDPA Problem;
// 	#ifdef DEBUG
//         cout << "making SDPA sparse format started " << endl;
// 	#endif
//         MakeSDPAform(sdpdata, Problem);
// 	#ifdef DEBUG
//        cout << "making SDPA sparse format finished " << endl;
// 	#endif
//         if(POP.param.SDPsolverOutFile.empty() == false){
//             fpresult = fopen(POP.param.SDPsolverOutFile.c_str(), "w");
//             if(fpresult == NULL){
//                 cout << " can not open file " << POP.param.SDPsolverOutFile << endl;
//             }
//             Problem.setResultFile(fpresult);
//         }
//         //cout << "setting file for writing SDPA result finished." << endl;
//         if(POP.param.SDPAOnScreen == 1){
//             Problem.setDisplay(stdout);
//         }else{
//             Problem.setDisplay(NULL);
//         }
//         //cout << "setting file for displaying SDPA result finished." << endl;
        
        
//         Problem.initializeSolve();
// 	    Problem.setParameterEpsilonStar(POP.param.SDPsolverEpsilon);
// 	    Problem.setParameterEpsilonDash(POP.param.SDPsolverEpsilon);
//         //cout << "solving SDP started." << endl;
//         Problem.solve();
//         //cout << "solving SDP finished " << endl;
//         //printSDPA(sdpdata);
//         cputime3 = (double)clock()/(double)CLOCKS_PER_SEC;
//         info.setandprintSDPInfo(sdpdata);
//         InputResultSDPA(Problem, info, POP.scalevalue[0], POP.Polysys.objConst);
//         double *y;
//         y = Problem.getResultXVec();
//         // Here is where y is created, so if we want to enforce any expectations equalities
//         // use expectMapping and do it here

// 	    int printSolSW = 1;
//         vector<double> xvect(POP.Polysys.dimVar);
//     	for(int i=0; i<POP.degOneTerms.size(); i++){
//         	if(POP.degOneTerms[i] == 0){
//         		printSolSW = 0;
// 			break;
// 		}
// 	}	
// 	if(printSolSW == 1){
// 	        for(int i=0; i<POP.Polysys.dimVar; i++){
//         	    xvect[i] = y[POP.degOneTerms[i]-1];
// 	        }

//         	// info.genAproximation(POP.Polysys.polynomial[0], POP.OriPolysys, POP.permmatrix, POP.bvect, xvect, fixedVar);
//             // Pass updated polysys instead of original one because we have changed it
//             info.genAproximation(POP.Polysys.polynomial[0], POP.Polysys, POP.permmatrix, POP.bvect, xvect, fixedVar);
//         	info.objValue = info.objValue*POP.scalevalue[0] + POP.Polysys.objConst;
// 		if(POP.Polysys.minOrmax == -1){
// 			info.objValue = -1 * info.objValue;
// 		}
//         }
//         cputime3 = (double)clock()/(double)CLOCKS_PER_SEC;
// 	    info.total = cputime3 - cputime1;
//         //cout << "infieasError = " << info.infeasError << endl;
// 	    //cout << "scaleError = " << info.scaleError << endl;
//         info.printSolutions(POP, POP.param.printFileName);
        
//         if(POP.param.printFileName.empty() == false && POP.param.printOnScreen == 1){
//             info.printSolutions(POP, "");
//         }
//         if(POP.param.SDPsolverOutFile.empty() == false){
//             fclose(fpresult);
//         }
//         Problem.terminate();
//     }else if(POP.param.SDPsolverSW == 1 && info.infeasibleSW == -1){
//         cout << "Not implemented, yet." << endl;
        
//     }else{
//         info.setandprintSDPInfo(sdpdata);
//         cputime3 = (double)clock()/(double)CLOCKS_PER_SEC;
//         info.total = cputime3 - cputime1;
//         info.printSolutions(POP, POP.param.printFileName);
//     }
//     return 0;
    
// }

