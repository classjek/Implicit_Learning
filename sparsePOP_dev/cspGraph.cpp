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
#include "cspGraph.h"

SymamdGraph::SymamdGraph(){
	ndim  = 0;
	cspnnz= 0;
	B     = NULL;
	q     = NULL;
	stats = NULL;
	perm  = NULL;
}

SymamdGraph::~SymamdGraph(){
	delete [] B;
	B = NULL;
	delete [] q;
	q = NULL;
	delete [] stats;
	stats = NULL;
	delete [] perm;
	perm   = NULL;
}

void SymamdGraph::makeGraph(vector<vector<double> > cspmat){
	stats = new int[COLAMD_STATS];
	ndim = cspmat[0].size();
	int nnz = 0;
	for(int i=0; i < ndim; i++){
		for(int j=i+1; j<ndim; j++){
			if(fabs(cspmat[j][i]) > EPS){
				nnz++;
			}
		}
	}
	cspnnz = nnz;
	//printCsp(cspmat);
	B = new int[cspnnz];
	q = new int[ndim+1];
	perm = new int[ndim+1];
	//int *stats = new int[COLAMD_STATS];
	int rownnz = 0;
	for(int i=0; i < ndim; i++){
		q[i] = rownnz;
		for(int j=i+1; j< ndim; j++){
			if(fabs(cspmat[j][i]) > EPS){
				B[rownnz] = j;
				rownnz++;
			}
		}
	}
	q[ndim] = cspnnz;
}
void SymamdGraph::disp(){
	printf("    B = ");
	for(int i=0; i<cspnnz; i++){
		printf("%2d ",B[i]);
	}      
	printf("\n");
	printf("    q = ");
	for(int i=0; i<ndim+1; i++){
		printf("%2d ",q[i]);
	}      
	printf("\n");
	printf(" perm = ");
	for(int i=0; i<ndim+1; i++){
		printf("%2d ",perm[i]);
	}      
	printf("\n");
}

MetisGraph::MetisGraph(){
	ndim = 0;
	numflag = 0;
	options= NULL;
	xadj   = NULL;
	adjncy = NULL;
	perm   = NULL;
	iperm  = NULL;
}

MetisGraph::~MetisGraph(){
	delete [] options;
	options = NULL;
	delete [] xadj;
	xadj   = NULL;
	delete [] adjncy;
	adjncy = NULL;
	delete [] perm;
	perm   = NULL;
	delete [] iperm;
	iperm  = NULL;
}

void MetisGraph::disp(){
	cout << " Information of Graph " << endl;
	cout << " ndim    = " << ndim << endl;
	cout << " numflag = " << numflag << endl;
	cout << " options " << endl;
	if(options){
		for(int i=0; i < 8; i++){
			cout << " options[" << i << "] = " << options[i] << endl;
		}
	}
	cout << " xadj " << endl;
	if(xadj){
		for(int i=0; i <= ndim; i++){
			cout << " xadj[" << i << "] = " << xadj[i] << endl;
		}
	}
	cout << " adjncy " << endl;
	if(adjncy){
		for(int i=0; i < xadj[ndim]; i++){
			cout << " adjncy[" << i << "] = " << adjncy[i] << endl;
		}
	}
	cout << " perm " << endl;
	if(perm){
		for(int i=0; i < ndim; i++){
			cout << " perm[" << i << "] = " << perm[i] << endl;
		}
	}
	cout << " iperm " << endl;
	if(iperm){
		for(int i=0; i < ndim; i++){
			cout << " iperm[" << i << "] = " << iperm[i] << endl;
		}
	}
	cout << endl;
}


void MetisGraph::makeGraph(vector<vector<double> > cspmat){

	options = new int[8];
	options[0] = 0;
	options[1] = 3;
	options[2] = 1;
	options[3] = 2;
	options[4] = 0;
	options[5] = 0;
	options[6] = 0;
	options[7] = 1;

	ndim = cspmat[0].size();
	int nnz  = 0;
	int pos  = 0;

	for(int i=0; i<ndim; i++){
		for(int j=i+1; j<ndim; j++){
			if(cspmat[i][j] != 0){
				nnz = nnz + 1;		
			}
		}
	}

	xadj   = new idxtype[ndim+1];
	adjncy = new idxtype[2*nnz];
	perm   = new idxtype[ndim];	
	iperm  = new idxtype[ndim];
	
	for(int i=0; i<ndim; i++){
		xadj[i] = pos;
		for(int j=0; j<i; j++){
			if(cspmat[i][j] != 0){
				adjncy[pos] = j;	
				pos = pos +1;
			}
		}
		for(int j=i+1; j<ndim; j++){
			if(cspmat[i][j] != 0){
				adjncy[pos] = j;	
				pos = pos +1;
			}
		}
	}
	xadj[ndim] = pos;
	//disp();
}

void MetisGraph::findPerm(){
	
	//METIS_NodeND(&ndim,xadj,adjncy,&numflag,options,perm,iperm);
	//disp();
}

void makeA(double * &A, vector<vector<double> > cspmat, vector<int> perm){
	int ndim = cspmat[0].size();
	int pos = 0;
	for(int i=0; i<ndim; i++){
		for(int j=i; j<ndim; j++){
			A[pos] = cspmat[perm[i]][perm[j]];
			pos = pos + 1;
		}
	}
} 
void makeExtmat(int nDim, double *A, class SparseMat & extmat){
	int pos = 0;
	int nnz = 0;
	
	vector<vector<int> > B;
	B.resize(nDim);
	for(int i=0; i<nDim; i++){
		B[i].resize(nDim,0);
	}

	for(int i=0;i<nDim; i++){
		for(int j=i; j<nDim; j++){
			if(fabs(A[pos]) >EPS){
				//cout << "A[" << pos << "]=" << A[pos] << endl;
				if(i != j){
					B[i][j] = 1;
					//B[j][i] = 1;
				}else{
					B[i][i] = 1;
				}
				nnz = nnz + 1;
			}
			pos = pos + 1;
		}
	}
	/*
	cout << "nnz = " << nnz << endl;
	cout << "pos = " << pos << endl;
	for(int i=0; i<nDim; i++){
		for(int j=0; j<nDim; j++){
			printf("%d ", B[i][j]);
		}
		printf("\n");
	}
	*/
	extmat.resizeIr(nnz, 0);
	extmat.resizeJc(nDim+1,0);	
	nnz = 0;
	for(int i=0;i<nDim; i++){
		extmat.jc[i] = nnz;
		for(int j=0; j<nDim; j++){
			if(B[i][j] == 1){
				extmat.ir[nnz] = j;
				nnz = nnz + 1;
			}
		}
	}
	extmat.jc[nDim] = nnz;
	/*
	for(int i=0; i<nnz;i++){
		cout << "extmat.ir[" << i << "] = " << extmat.ir[i] << endl;
	}
	for(int i=0; i<nDim+1;i++){
		cout << "extmat.jc[" << i << "] = " << extmat.jc[i] << endl;
	}
	*/
}
/*
void makeExtmat(double *A, vector<vector<int> > & extmat){
	int pos = 0;
	int ndim = extmat[0].size();
	for(int i=0; i<ndim; i++){
		//printf("ndim = %2d\n", extmat[i].size());
		for(int j=i; j<ndim; j++){
			if(fabs(A[pos]) > EPS){
				//printf("(i,j) = (%2d, %2d)\n", i,j);
				extmat[i][j] = 1;
				extmat[j][i] = 1;
			}
			pos = pos + 1;
		}
	}
}
*/
void genCSP(class polysystem Polysys, vector<vector<double> > & cspmat){
	
	cspmat.resize(Polysys.dimVar);
	for(int i=0; i<Polysys.dimVar; i++){
		cspmat[i].resize(Polysys.dimVar,0);
		//rv = (double)(rand()%100000)/(100000);
		//while(fabs(rv) < EPS){
		//	rv = (double)(rand()%100000)/(100000);
		//}
		//rv = (double)(1+rand()%100000)/(100000);
		//cspmat[i][i] = cspmat[i][i] + rv; 
		cspmat[i][i] = (double)2*6*Polysys.dimVar;
	}	
	
	list<mono>::iterator it;
	vector<poly> polyVec = Polysys.polynomial;
	vector<int> comb,sumSupp;
	int deg;
	for(int i=0; i < Polysys.numSys; i++){
		if(i == 0){
			for(it = polyVec[0].monoList.begin(); it != polyVec[0].monoList.end();++it){
				//deg = accumulate((*it).Supp.begin(),(*it).Supp.end(),0);
				deg = accumulate((*it).supVal.begin(),(*it).supVal.end(),0);
				if(deg > 1){
					(*it).getSuppComb(comb);
					intVecUnique(comb);
					assignCmat(cspmat,comb);
				}
				/*
				cout << "i = " << 0 << endl;
				for(int j=0; j < comb.size(); j++){
					cout << "comb[" << j << "]= " << comb[j] << endl;
				}
				*/
				comb.clear();
			}
		}else{
			polyVec[i].sumSupports(sumSupp);
			for(int j = 0; j < sumSupp.size(); j++){
				//if(sumSupp[j] != 0){
					int v = sumSupp[j];
					comb.push_back(v);
				//}
			}
			intVecUnique(comb);
			assignCmat(cspmat,comb);
			/*
			cout << "i = " << i << endl;
			for(int j=0; j < comb.size(); j++){
				cout << "comb[" << j << "]= " << comb[j] << endl;
			}
			*/
			comb.clear();
		}
	}
	for(int i=0; i<Polysys.dimVar; i++){
		//cspmat[i][i] = 6*Polysys.dimVar*(double)(1+rand()%100000)/100000;
		cspmat[i][i] = cspmat[i][i] + 5*Polysys.dimVar;
	}	
}

void assignCmat(vector<vector<double> > & cspmat, vector<int> comb){
	vector<int>::iterator it,it1;
	double rv;
	for(it = comb.begin(); it != comb.end(); ++it){
		for(it1 = it+1; it1 != comb.end(); ++it1){
			//rv = (double)(rand()%100000)/(100000);
			//while(fabs(rv) < EPS){
			//	rv = (double)(rand()%100000)/(100000);
			//}
			rv = (double)(1+rand()%100000)/(100000);
			cspmat[(*it)][(*it1)] = rv;
			cspmat[(*it1)][(*it)] = rv;
		}
	}
}

void intVecUnique(vector<int> & ivec){
	set<int> ret;
	for(vector<int>::iterator it = ivec.begin(); it != ivec.end(); ++it){
		ret.insert(*it);
	}
	ivec.clear();
	for(set<int>::iterator it = ret.begin();it != ret.end(); ++it){
		ivec.push_back(*it);
	}
}

