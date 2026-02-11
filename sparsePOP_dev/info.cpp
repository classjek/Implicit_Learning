/* -------------------------------------------------------------

This file is a component of SparsePOP
Copyright (C) 2007 SparsePOP Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */
#include "info.h"

Info::Info(){
	sizeArow    = 0;
	sizeAcol    = 0;
	Annz        = 0;
	kf          = 0;
	kl          = 0;
	maxks       = 0;
	aveks       = 0;
	infeasibleSW= 0;
	reduceAMatSW= 0;
	SDPobjValue = 0;
	objValue    = 0;
	infeasError = 0;
	scaleError  = 0;
	conversion  = 0;
	SDPAtime    = 0;
	total       = 0;
	phase       = new char[15];
}
Info::~Info(){
	ks.clear();
	xvect.clear();
	delete [] phase;
}

void Info::printSolutions(class s3r SDPr, string fname){
	if(fname.empty() == true && SDPr.param.printOnScreen == 1){
		printf("\n\n## Computational Results by sparsePOP with SDPA ##\n");
		//printf("## Printed by printSolustions ##\n");
		printf("# Problem File Name   = %s\n",SDPr.problemName.c_str());	
		printf("# parameters:\n");
		printf("  relaxOrder          = %d\n", SDPr.param.relax_Order);
		//if(SDPr.param.sparseSW == 0){
			printf("  sparseSW            = %d\n", SDPr.param.sparseSW);
		//}
		if(SDPr.param.multiCliquesFactor != 1){
			printf("  multiCliquesFactor  = %d\n", (int)(SDPr.param.multiCliquesFactor));
		}
		if(SDPr.param.scalingSW == 0){
			printf("  scalingSW           = %d\n", SDPr.param.scalingSW);
		}
		if(SDPr.param.boundSW != 2){
			printf("  boundSW             = %d\n", SDPr.param.boundSW);
		}
		if(fabs(SDPr.param.eqTolerance) > EPS*0.01){
			printf("  eqTolerance         = %6.2e\n", SDPr.param.eqTolerance);
		}
		if(fabs(SDPr.param.perturbation) > EPS*0.01){
			printf("  perturbation        = %6.2e\n", SDPr.param.perturbation);
		}
		if(SDPr.param.reduceMomentMatSW != 1){
			printf("  reduceMomentMatSW   = %d\n", SDPr.param.reduceMomentMatSW);
		}
		if(SDPr.param.complementaritySW == 1){
			printf("  complementaritySW   = %d\n", SDPr.param.complementaritySW);
    		}
		if(SDPr.param.SquareOneSW == 0){
			printf("  SquareOneSW         = %d\n", SDPr.param.SquareOneSW);
    		}
		if(SDPr.param.binarySW == 0){
			printf("  binarySW            = %d\n", SDPr.param.binarySW);
    		}
		if(SDPr.param.SDPsolverOutFile.empty() == false){
			printf("  SDPsolverOutFile    = %s\n",SDPr.param.SDPsolverOutFile.c_str());	
    		}
		if(SDPr.param.detailedInfFile.empty() == false){
			printf("  detailedInfFile     = %s\n",SDPr.param.detailedInfFile.c_str());	
    		}
		if(SDPr.param.sdpaDataFile.empty() == false){
			printf("  sdpaDataFile        = %s\n",SDPr.param.sdpaDataFile.c_str());	
    		}
	    	if(SDPr.param.printFileName.empty() == false){
			printf("  printFileName       = %s\n",SDPr.param.printFileName.c_str());	
 	   	}else if(SDPr.param.printOnScreen == 0 && SDPr.param.printFileName.empty() == true){
    	    		printf("  printFileName       = 0\n");
    		}
		if(SDPr.param.Method.empty() == false){
			printf("  Method              = %s\n",SDPr.param.Method.c_str());	
    		}

		if(infeasibleSW >= -1 && infeasibleSW <= 1){
			if(SDPr.param.SDPsolverSW == 1 && infeasibleSW == 0){
				printf("# SDP solved by SDPA \n");
			}else if(SDPr.param.SDPsolverSW == 0 && infeasibleSW == 0){
				printf("# Estimated size of SDP to be solved: \n");
			}else if( reduceAMatSW <= 1){
				printf("# SDP relaxation problem: \n");
			}
			if(reduceAMatSW <= 1){
				printf("  size of A           = [%d, %d]\n", sizeArow, sizeAcol);
				printf("  no of nonzeros in A = %d\n", Annz);
				printf("  no of LP variables  = %d\n", kl);  
				printf("  no of FR variables  = %d\n",kf);  
				printf("  no of SDP blocks    = %lu\n", ks.size());  
				printf("  max size SDP block  = %d\n", maxks);  
				printf("  ave.size SDP block  = %3.2f\n", aveks); 
			} 
		}	
		if((SDPr.param.SDPsolverSW == 1 && infeasibleSW == 0)||infeasibleSW == -1){
			printf("# SDPA information \n");
			// string::iterator it = remove_if(sdpaPhase.begin(), sdpaPhase.end(),bind2nd(equal_to<char>(),' '));
			string::iterator it = remove_if(sdpaPhase.begin(), sdpaPhase.end(), [](char c){ return c == ' '; });
			sdpaPhase = string(sdpaPhase.begin(), it);
			printf("          phase value = %s\n",sdpaPhase.c_str());	
			
			/* SDPA information */
			printf("            Iteration = %d\n",iter);
			printf("                   mu = %+8.5e\n", mu);
			printf("         relative gap = %+8.5e\n", relativeGap);
			printf("                  gap = %+8.5e\n", gap);
			printf("               digits = %+8.5e\n", digit);
			//printf("         objValPrimal = %+8.5e\n",PrimalObj);
			//printf("         objValDual   = %+8.5e\n", DualObj);
			printf("         Primal Error = %+8.5e\n", PrimalError);
			printf("           Dual Error = %+8.5e\n", DualError);

			if(sdpaPhase == "pdINF"){
				printf(" At lease one of the primal and the dual SDP is expected to be infeasible.\n");
			}else if(sdpaPhase == "pINF"||sdpaPhase == "pINF_dFEAS"){
				printf(" Primal SDP is infeasible.\n");
				printf(" POP is probably infeasible.\n");
			}else if(sdpaPhase == "dINF"||sdpaPhase == "pFEAS_dINF"){
				printf(" Dual SDP is infeasible.\n");
				printf(" No finite lower bound has been computed for POP objective function.\n");
			}else if(sdpaPhase == "noINFO"){
				printf(" The iteration has exceeded the maxIterationx and/or stopped with\n");
				printf(" no information on the primal and dual feasibility.\n");
			}else if(sdpaPhase  == "pUNBD"){
				printf(" Primal SDP may be unbounded.\n");
				printf(" No finite lower bound has been computed for POP objective function.\n");
			}else if(sdpaPhase == "dUNBD"){
				printf(" Dual SDP may be unbounded.\n");
				printf(" POP is probably infeasible.\n");
			}
		}
		
		if(infeasibleSW == -1|| sdpaPhase == "pdOPT" || sdpaPhase == "pFEAS" || sdpaPhase == "dFEAS" || sdpaPhase == "pdFEAS"){/* I will use SDPA information */
			printf("# Approximate optimal value information:\n");
			if(infeasibleSW == 0){
				printf("  SDPobjValue         = %+13.7e\n", SDPobjValue);
			}
			if(xvect.empty() == false){
				printf("  POP.objValue        = %+13.7e\n", objValue);
				if(infeasibleSW == 0){
					double ma = 1;
					if(ma < fabs(objValue)){
						ma = fabs(objValue);
					}
					double relobj = fabs(objValue - SDPobjValue)/ma;
					printf("  relatvie obj error  = %+8.3e\n", relobj);
				}
				printf("  POP.absError        = %+8.3e\n", infeasError);
				printf("  POP.scaledError     = %+8.3e\n", scaleError);
			}
		}else if(infeasibleSW == 1){
			printf("# Dual SDP is infeasible!\n");
			printf("  No finite lower bound has been computed for the POP objective function.\n");
		}else if(infeasibleSW == 2){
			printf("# POP is infeasible!\n");
		}else if(infeasibleSW == -2){
			printf("# Approximate optimal value information:\n");
			printf("  POP.objValue        = %+13.7e\n", objValue);
			printf("  POP.absError        = %+8.3e\n", infeasError);
			printf("  POP.scaledError     = %+8.3e\n", scaleError);
		}
		/**/	
		printf("# cpu time: \n");
		printf("  cpuTime.conversion  = %8.2f\n", conversion);
		if(infeasibleSW == 0 && SDPr.param.SDPsolverSW == 1){
			printf("  cpuTime.SDPA        = %8.2f\n", SDPAtime);
		}
		printf("  cpuTime.total       = %8.2f\n", total);
		if(SDPr.param.SDPsolverSW == 1 && xvect.empty() == false && SDPr.param.printLevel1 >= 2){
			if(sdpaPhase == "pdOPT" || sdpaPhase == "pFEAS" || sdpaPhase == "dFEAS" || sdpaPhase == "pdFEAS"|| infeasibleSW == -2){
				printf("# Approximate optimal solution information:\n");
				for(int i=0; i<xvect.size(); i++){
					if(i % 5 == 0){
						printf(" ");
					}	
					if((i+1) % 5 == 0){
						printf(" %4d: %+13.7e\n",i+1,xvect[i]);
					}else{
						printf(" %4d: %+13.7e",i+1,xvect[i]);
					}
				}
				printf("\n");
			}
		}else if(SDPr.param.SDPsolverSW == 1 && xvect.empty() == true && SDPr.param.printLevel1 >= 2){
			cout << "Oh no are we here!? Yikes" << endl;
			printf("\n## SparsePOP can not return approximate\n");
			printf("## optimal solution information for this problem.\n");
			printf("## If you need this solutino information, \n");
			printf("## please solve with param.reduceMomentMatSW = 0.\n");
		}
	}else if(fname.empty() == false){
		
		FILE *fp;
		fp = fopen(fname.c_str(),"a+");
		if(fp == NULL){
			printf("file open error\n");
			exit(EXIT_FAILURE);
		}
		fprintf(fp,"## Computational Results by sparsePOP with SDPA ##\n");
		//fprintf(fp,"## Printed by printSolustions ##\n");
		fprintf(fp,"# Problem File Name   = %s\n",SDPr.problemName.c_str());	
		fprintf(fp,"# parameters:\n");
		fprintf(fp,"  relaxOrder          = %d\n", SDPr.param.relax_Order);
		//if(SDPr.param.sparseSW == 0){
			fprintf(fp,"  sparseSW            = %d\n", SDPr.param.sparseSW);
		//}
		if(SDPr.param.multiCliquesFactor != 1){
			fprintf(fp,"  multiCliquesFactor  = %d\n", (int)(SDPr.param.multiCliquesFactor));
		}
		if(SDPr.param.scalingSW == 0){
			fprintf(fp,"  scalingSW           = %d\n", SDPr.param.scalingSW);
		}
		if(SDPr.param.boundSW != 2){
			fprintf(fp,"  boundSW             = %d\n", SDPr.param.boundSW);
		}
		if(fabs(SDPr.param.eqTolerance) > EPS*0.01){
			fprintf(fp,"  eqTolerance         = %6.2e\n", SDPr.param.eqTolerance);
		}
		if(fabs(SDPr.param.perturbation) > EPS*0.01){
			fprintf(fp,"  perturbation        = %6.2e\n", SDPr.param.perturbation);
		}
		if(SDPr.param.reduceMomentMatSW == 0){
			fprintf(fp,"  reduceMomentMatSW   = %d\n", SDPr.param.reduceMomentMatSW);
		}
		if(SDPr.param.complementaritySW == 1){
			fprintf(fp,"  complementaritySW   = %d\n", SDPr.param.complementaritySW);
    		}
		if(SDPr.param.SquareOneSW == 0){
			fprintf(fp,"  SquareOneSW         = %d\n", SDPr.param.SquareOneSW);
    		}
		if(SDPr.param.binarySW == 0){
			fprintf(fp,"  binarySW            = %d\n", SDPr.param.binarySW);
    		}
		if(SDPr.param.SDPsolverOutFile.empty() == false){
			fprintf(fp,"  SDPsolverOutFile    = %s\n",SDPr.param.SDPsolverOutFile.c_str());	
    		}
		if(SDPr.param.detailedInfFile.empty() == false){
			fprintf(fp,"  detailedInfFile     = %s\n",SDPr.param.detailedInfFile.c_str());	
    		}
		if(SDPr.param.sdpaDataFile.empty() == false){
			fprintf(fp,"  sdpaDataFile        = %s\n",SDPr.param.sdpaDataFile.c_str());	
    		}
		if(SDPr.param.detailedInfFile.empty() == false){
			fprintf(fp,"  detailedInfFile     = %s\n",SDPr.param.detailedInfFile.c_str());	
    		}
		if(SDPr.param.sdpaDataFile.empty() == false){
			fprintf(fp,"  sdpaDataFile        = %s\n",SDPr.param.sdpaDataFile.c_str());	
		}
	    	if(SDPr.param.printFileName.empty() == false){
			fprintf(fp,"  printFileName       = %s\n",SDPr.param.printFileName.c_str());
 	   	}else if(SDPr.param.printOnScreen == 0 && SDPr.param.printFileName.empty() == true){
    	    		fprintf(fp,"  printFileName       = 0\n");
    		}
		if(SDPr.param.Method.empty() == false){
			fprintf(fp,"  Method              = %s\n",SDPr.param.Method.c_str());	
    		}

		if(infeasibleSW >= -1 && infeasibleSW <= 1){
			if(SDPr.param.SDPsolverSW == 1 && infeasibleSW == 0){
				fprintf(fp,"# SDP solved by SDPA \n");
			}else if(SDPr.param.SDPsolverSW == 0 && infeasibleSW == 0){
				fprintf(fp,"# Estimated size of SDP to be solved: \n");
			}else if( reduceAMatSW <= 1){
				fprintf(fp,"# SDP relaxation problem: \n");
			}
			if(reduceAMatSW <= 1){
				fprintf(fp,"  size of A           = [%d, %d]\n", sizeArow, sizeAcol);
				fprintf(fp,"  no of nonzeros in A = %d\n", Annz);
				fprintf(fp,"  no of LP variables  = %d\n", kl);  
				fprintf(fp,"  no of FR variables  = %d\n",kf);  
				fprintf(fp,"  no of SDP blocks    = %lu\n", ks.size());  
				fprintf(fp,"  max size SDP block  = %d\n", maxks);  
				fprintf(fp,"  ave.size SDP block  = %3.2f\n", aveks); 
			} 
		}	
		if((SDPr.param.SDPsolverSW == 1 && infeasibleSW == 0)||infeasibleSW == -1){
			fprintf(fp,"# SDPA information \n");
			// string::iterator it = remove_if(sdpaPhase.begin(), sdpaPhase.end(),bind2nd(equal_to<char>(),' '));
			string::iterator it = remove_if(sdpaPhase.begin(), sdpaPhase.end(), [](char c){ return c == ' '; });
			sdpaPhase = string(sdpaPhase.begin(), it);
			fprintf(fp,"          phase value = %s\n",sdpaPhase.c_str());
			/* SDPA information */
			fprintf(fp,"            Iteration = %d\n",iter);
			fprintf(fp,"                   mu = %+8.5e\n", mu);
			fprintf(fp,"         relative gap = %+8.5e\n", relativeGap);
			fprintf(fp,"                  gap = %+8.5e\n", gap);
			fprintf(fp,"               digits = %+8.5e\n", digit);
			//printf("         objValPrimal = %+8.5e\n",PrimalObj);
			//printf("         objValDual   = %+8.5e\n", DualObj);
			fprintf(fp,"         Primal Error = %+8.5e\n", PrimalError);
			fprintf(fp,"         Dual Error   = %+8.5e\n", DualError);

			if(sdpaPhase == "pdINF"){
				fprintf(fp," At lease one of the primal and the dual SDP is expected to be infeasible.\n");
			}else if(sdpaPhase == "pINF"||sdpaPhase == "pINF_dFEAS"){
				fprintf(fp," Primal SDP is infeasible.\n");
				fprintf(fp," POP is probably infeasible.\n");
			}else if(sdpaPhase == "dINF"||sdpaPhase == "pFEAS_dINF"){
				fprintf(fp," Dual SDP is infeasible.\n");
				fprintf(fp," No finite lower bound has been computed for POP objective function.\n");
			}else if(sdpaPhase == "noINFO"){
				fprintf(fp," The iteration has exceeded the maxIterationx and stopped with\n");
				fprintf(fp," no information on the primal and dual feasibility.\n");
			}else if(sdpaPhase  == "pUNBD"){
				fprintf(fp," Primal SDP may be unbounded.\n");
				fprintf(fp," No finite lower bound has been computed for POP objective function.\n");
			}else if(sdpaPhase == "dUNBD"){
				fprintf(fp," Dual SDP may be unbounded.\n");
				fprintf(fp," POP is probably infeasible.\n");
			}
		}
		
		if(infeasibleSW == -1|| sdpaPhase == "pdOPT" || sdpaPhase == "pFEAS" || sdpaPhase == "dFEAS" || sdpaPhase == "pdFEAS"){/* I will use SDPA information */
			fprintf(fp,"# Approximate optimal value information:\n");
			if(infeasibleSW == 0){
				fprintf(fp,"  SDPobjValue         = %+13.7e\n", SDPobjValue);
			}
			if(xvect.empty() == false){
				fprintf(fp,"  POP.objValue        = %+13.7e\n", objValue);
				if(infeasibleSW == 0){
					double ma = 1;
					if(ma < fabs(objValue)){
						ma = fabs(objValue);
					}
					double relobj = fabs(objValue - SDPobjValue)/ma;
					fprintf(fp,"  relatvie obj error  = %+8.3e\n", relobj);
				}
				fprintf(fp,"  POP.absError        = %+8.3e\n", infeasError);
				fprintf(fp,"  POP.scaledError     = %+8.3e\n", scaleError);
			}
		}else if(infeasibleSW == 1){
			fprintf(fp,"# Dual SDP is infeasible!\n");
			fprintf(fp,"  No finite lower bound has been computed for the POP objective function.\n");
		}else if(infeasibleSW == 2){
			fprintf(fp,"# POP is infeasible!\n");
		}else if(infeasibleSW == -2){
			fprintf(fp,"# Approximate optimal value information:\n");
			fprintf(fp,"  POP.objValue        = %+13.7e\n", objValue);
			fprintf(fp,"  POP.absError        = %+8.3e\n", infeasError);
			fprintf(fp,"  POP.scaledError     = %+8.3e\n", scaleError);
		}
		/**/	
		fprintf(fp,"# cpu time: \n");
		fprintf(fp,"  cpuTime.conversion  = %8.2f\n", conversion);
		if(infeasibleSW == 0 && SDPr.param.SDPsolverSW == 1){
			fprintf(fp,"  cpuTime.SDPA        = %8.2f\n", SDPAtime);
		}
		fprintf(fp,"  cpuTime.total       = %8.2f\n", total);
		if(SDPr.param.SDPsolverSW == 1 && xvect.empty() == false && SDPr.param.printLevel2 >= 2){
			if(sdpaPhase == "pdOPT" || sdpaPhase == "pFEAS" || sdpaPhase == "dFEAS" || sdpaPhase == "pdFEAS"|| infeasibleSW == -2){
				fprintf(fp,"# Approximate optimal solution information:\n");
				for(int i=0; i<xvect.size(); i++){
					if(i % 5 == 0){
						fprintf(fp," ");
					}	
					if((i+1) % 5 == 0){
						fprintf(fp," %4d: %+13.7e\n",i+1,xvect[i]);
					}else{
						fprintf(fp," %4d: %+13.7e",i+1,xvect[i]);
					}
				}
				fprintf(fp,"\n");
			}
		}else if(SDPr.param.SDPsolverSW == 1 && xvect.empty() == true && SDPr.param.printLevel2 >= 2){
			fprintf(fp, "\n## SparsePOP can not return approximate\n");
			fprintf(fp, "## optimal solution information for this problem.\n");
			fprintf(fp, "## If you need this solutino information, \n");
			fprintf(fp, "## please solve with param.reduceMomentMatSW = 0.\n");
		}
		fprintf(fp,"\n");
		fclose(fp);
	}

}

void Info::CheckFeasibility(class polysystem Polysys, double & infeasError, double & scaleError){
	vector<double> fValue;
	vector<double> maxAbs; 
	/*
	Polysys.polynomial[0].evalPoly(xvect,fValue,maxAbs);
	objValue = fValue[0];
	Polysys.polynomial[0].writePolyData();
	*/

	int size;
	double maxA;
	for(int i=1; i<Polysys.numSys; i++){
		fValue.clear();
		maxAbs.clear();
		Polysys.polynomial[i].evalPoly(xvect,fValue,maxAbs);
		//Polysys.polynomial[i].writePolyData();
		if(Polysys.polynomial[i].typeCone == 1){
			size = Polysys.polynomial[i].sizeCone;
			for(int i=0; i<size; i++){
				maxA = 1;
				if(maxA < maxAbs[i]){
					maxA = maxAbs[i];
				}
				if(infeasError > fValue[i]){
					infeasError = fValue[i];
				}
				if(scaleError > fValue[i]/maxA){
					scaleError = fValue[i]/maxA;
				}
			}
		}else if(Polysys.polynomial[i].typeCone == -1){
			size = Polysys.polynomial[i].sizeCone;
			for(int i=0; i<size; i++){
				maxA = 1;
				if(maxA < maxAbs[i]){
					maxA = maxAbs[i];
				}
				if(infeasError > -fabs(fValue[i])){
					infeasError = -fabs(fValue[i]);
				}
				if(scaleError > -fabs(fValue[i]/maxA)){
					scaleError = -fabs(fValue[i]/maxA);
				}
			}

		}else{
			cout << " ## Error: Not defined " << endl;
			exit(EXIT_FAILURE);
		}	
		/*	
		cout << "Polynomial " << i << endl;
		cout << "fValue = " << fValue[0] << endl;
		cout << "maxAbs = " << maxAbs[0] << endl;
		cout << "infeasError = " << infeasError << endl;
		cout << "scaleError  = " << scaleError << endl;
		*/
	}

	/* check bounds */
	double ltmp  = 0;
	double mltmp = 0;
	double utmp  = 0;
	double mutmp = 0;
	for(int i=0; i<Polysys.dimVar; i++){
		if(Polysys.bounds.lbd(i) > MIN){
			ltmp = xvect[i] - Polysys.bounds.lbd(i);
			if(fabs(Polysys.bounds.lbd(i)) > fabs(xvect[i])){ 
				mltmp = fabs(Polysys.bounds.lbd(i)); 
			}else{
				mltmp = fabs(xvect[i]);
			}	
			if(mltmp < 1){
				mltmp = 1;
			}
			if(infeasError > ltmp){
				infeasError = ltmp;
			}
			if(scaleError > ltmp/mltmp){
				scaleError = ltmp/mltmp;
			}
		}
		if(Polysys.bounds.ubd(i) < MAX){
			utmp = Polysys.bounds.ubd(i) - xvect[i];
			if(fabs(Polysys.bounds.ubd(i)) > fabs(xvect[i])){ 
				mutmp = fabs(Polysys.bounds.ubd(i)); 
			}else{
				mutmp = fabs(xvect[i]);
			}
			if(mutmp < 1){
				mutmp = 1;
			}
			if(infeasError > utmp){
				infeasError = utmp;
			}
			if(scaleError > utmp/mutmp){
				scaleError = utmp/mutmp;
			}
		}
	}
}

void Info::nDimZero(class polysystem Polysys, vector<vector<double> > fixedVar){
	reduceAMatSW = 2;
	xvect.resize(fixedVar[1].size());
	for(int i=0; i<fixedVar[1].size(); i++){
		xvect[i] = fixedVar[1][i];
		//cout << xvect[i] << " ";
	}
	//cout << endl;
	double infeasError = 0;
	double scaleError = 0;
	CheckFeasibility(Polysys, infeasError, scaleError);
	if(infeasError < -EPS || scaleError < -EPS){	
		infeasibleSW = 2;
		cout << "## POP is infeasible." << endl;
	}else{
		infeasibleSW = -2;
		cout << "## POP has a unique feasible solution." << endl;
	}
}

// if we have extended fixedVar to include added variables 
// then hopefully this should go smoothly
void Info::genAproximation(class poly objPoly, class polysystem Polysys, vector<double> permmatrix, vector<double> bvect, vector<double> vect, vector<vector<double> > fixedVar){
	vector<double> fValue;
	vector<double> maxAbs; 

	objPoly.evalPoly(vect,fValue,maxAbs);
	objValue = fValue[0];

	// Project variables back //
	vector<int> idx; // collect indices of variables that were not eliminated
	xvect.resize(fixedVar[0].size()); 
	// fill in eliminated variables and collect remaining indices
	for(int i=0; i<fixedVar[0].size(); i++){ 
		if(fixedVar[0][i] == 1){
			xvect[i] = fixedVar[1][i];
		}else{
			idx.push_back(i);
		}
	}

	// Map reduced solutions back to original variables
	int ndim = idx.size();
	for(int i=0; i<ndim; i++){	
		xvect[idx[i]] = vect[i]*permmatrix[i]  + bvect[i]; // apply affine transformation (scaling/shifting)
	}
	// if variable i was eliminated during preprocessing(set to constant)
	// then xvect[i] holds the constant value that the variable was fixed to
	// if variable i survived preprocessing, xvect[i] holds the value of the varaible 
	// as computed from the solution
	
	// Debug 
	cout << "xvect(after scale) = [";
	for (size_t i = 0; i < xvect.size(); ++i) {
    	cout << xvect[i];
    	if (i + 1 < xvect.size()) cout << ", ";
	}
	cout << "]" << endl;

	infeasError = 0;
	scaleError = 0;
	CheckFeasibility(Polysys, infeasError, scaleError);
	cout << "infeasError = " << infeasError << endl;  // debug
	cout << "scaleError  = " << scaleError << endl;  // debug
} 

void Info::setandprintSDPInfo(class mysdp sdpdata){
	ks.empty();
	kl       = 0;
	kf       = 0;
	sizeArow = sdpdata.mDim;
	sizeAcol = 0;
	Annz     = 0;	
	maxks    = 0;
	aveks    = 0;
	SDPobjValue = 0;

	int kssize = 0;
	int tmp;
	for(int i=1; i<sdpdata.nBlocks+1;i++){
		tmp = sdpdata.bLOCKsTruct[i];
		if(sdpdata.block_info[2][i] == 2 || sdpdata.block_info[2][i] == 3){
			tmp = tmp*tmp;
			kssize = kssize + 1;
		}else if(sdpdata.block_info[2][i] == 1){
			tmp = -tmp;
			kl = kl + tmp;
		}else if(sdpdata.block_info[2][i] == -1){
			tmp = -tmp;
			tmp = (tmp >> 1);
			kf = kf + tmp;
		}
		sizeAcol = sizeAcol + tmp;
	}
	int tmpAnnz = 0;
	for(int i=1; i<sdpdata.nBlocks+1;i++){
		tmp = sdpdata.block_info[0][i] + sdpdata.block_info[1][i];
		for(int j=sdpdata.block_info[0][i];j<tmp;j++){
			if(sdpdata.ele.sup.pnz[0][j] != 0){/* Remove nnz of C */
				if(sdpdata.block_info[2][i] == -1){/* typeCone = -1 */
					tmpAnnz = tmpAnnz + 1;
				}else{
					if(sdpdata.ele.bij[1][j] == sdpdata.ele.bij[2][j]){
						Annz = Annz + 1;		
					}else{
						Annz = Annz + 2;
					}
				}
			}
		}
		if(sdpdata.block_info[2][i] == 2 || sdpdata.block_info[2][i] == 3){
			ks.push_back(sdpdata.bLOCKsTruct[i]);
		}
	}
	tmpAnnz = (tmpAnnz >> 1);
	Annz = Annz + tmpAnnz;
	vector<int>::iterator it;
	maxks = 0;
	aveks = 0;
	for(it = ks.begin(); it != ks.end(); ++it){
		if(maxks < (*it)){
			maxks = (*it);
		}
		aveks = aveks + (*it);
	}
	if(kssize >0){
		aveks = aveks/kssize;
	}	
	/*
	cout << "size(A) = [" << sizeArow << ", " << sizeAcol << "];" << endl;
	cout << "nnz(A)  = " << Annz << endl;
	cout << "K.l     = " << kl << endl;
	cout << "K.f     = " << kf << endl;
	cout << "K.s     = ";	
	for(it = ks.begin(); it != ks.end(); ++it){
		cout << *it << " ";
	}
	cout << endl;
	cout << "maxks   = " << maxks << endl;
	cout << "aveks   = " << aveks << endl;
	*/
}

