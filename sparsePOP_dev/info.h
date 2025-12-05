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
#ifndef __info_h__
#define __info_h__
#include "global.h"
#include "conversion.h"
#include "polynomials.h"
#include "mysdp.h"

class Info{
	/* information for case of SeDuMi */
	public:
		int sizeArow;
		int sizeAcol;
		int Annz;
		int kf;
		int kl;
		vector<int> ks;
		int maxks;
		double aveks;
	
		int infeasibleSW;
		/*
		infeasibleSW = 2; POP is found infeasible before applying SDPA.
		infeasibleSW = 1; Primal SDP is infeasible before applying SDPA.
		infeasibleSW = 0; Alppy the sparse SDP relaxation.
		infeasibleSW = -1; SDP is sovled by preprocessor of sparsePOP.
		infeasibleSW = -2; POP has a unique feasible soulution.
		*/	
		int reduceAMatSW;
		/*
		reduceAMatSW = 2; all the variables of POP are fixed.
		reduceAMatSW = 1; some equality constraints of the SDPA format dual SDP are lienarly dependant.
		reduceAMatSW = 0; no variable of POP is fixed.
		*/
		/* SDPA information */
		double PrimalObj;
		double DualObj;
		double PrimalError;
		double DualError;
		char*  phase;		
		string sdpaPhase;
		double relativeGap;
		double gap;
		double mu;
		int iter;
		double digit;
			
		/* cpu times */
		double conversion;
		double SDPAtime;
		double total;
	
	/* results obtained by SDP relaxation */
		double SDPobjValue;
		double objValue;
		double infeasError;
		double scaleError;
		vector<double> xvect;
	
		Info();
		~Info();	
	
		void setandprintSDPInfo(class mysdp sdpdata);
		void genAproximation(class poly objPoly, class polysystem Polysys, vector<double> permmatrix, vector<double> bvect, vector<double> xvect, vector<vector<double> > fixedVar); 
		void printSolutions(class s3r, string fname);
		void nDimZero(class polysystem Polysys, vector<vector<double> > fixedVar);
	
		void CheckFeasibility(class polysystem Polysys, double & infeasibleError, double & scaleError);
};

#endif // __info_h__
