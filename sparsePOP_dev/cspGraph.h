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
#ifndef __cspGraph_h__
#define __cspGraph_h__
#include "global.h"
#include "polynomials.h"
#include "conversion.h"

void makeA(double * &, vector<vector<double> >, vector<int>);
void makeExtmat(int, double *, class SparseMat &);
void makeExtmat(double *A, vector<vector<int> > & extmat);
void genCSP(class polysystem, vector<vector<double> > &);
void assignCmat(vector<vector<double> > &, vector<int>);
void intVecUnique(vector<int> &);

class MetisGraph{
	public:
		int ndim;
		int numflag;
		int  *options;
		idxtype *xadj;
		idxtype *adjncy;
		idxtype *perm;
		idxtype *iperm;

	MetisGraph();
	~MetisGraph();
	void makeGraph(vector<vector<double> > cspmat);
	void findPerm();
	void disp();

};

class SymamdGraph{
	public:
		int ndim;
		int cspnnz;
		int *B;
		int *q;
		int *stats;
		int *perm;

	SymamdGraph();
	~SymamdGraph();
	void makeGraph(vector<vector<double> > cspmat);
	void disp();
};

#endif // __cspGraph_h__
