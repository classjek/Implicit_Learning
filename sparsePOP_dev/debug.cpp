#include "debug.h"

void printDegOne(vector<int> vec){
	cout << "degOne = ";
	vector<int>::iterator it;
	for(it = vec.begin();it != vec.end(); ++it){
		cout << *it << " ";
	} 
	cout << endl;
}

void printSDPA(class mysdp sdpdata){
	cout << " --- sdpdata --- " << endl;
    cout << "mDim    = " << sdpdata.mDim << endl;
    cout << "nBlocks = " << sdpdata.nBlocks << endl;
    cout << "bLOCKsTruct " << endl;
	for(int i=0; i<sdpdata.bLOCKsTruct.size(); i++){
        cout << sdpdata.bLOCKsTruct[i] << " ";
    }
    cout << endl;
    cout << "block_info" << endl;
    for(int i=0; i < 3; i++){
        cout << "block_info[" << i << "] = ";
        for(int j=0; j < sdpdata.block_info[i].size(); j++){
            cout << sdpdata.block_info[i][j] << " ";
        }
        cout << endl;
    }
    cout << " mat_info" << endl;
    cout << " mat_info.bij" << endl;
    for(int i=0; i < 3; i++){
        cout << "mat_info.bij[" << i << "] = ";
        for(int j=0; j < sdpdata.ele.bij[i].size(); j++){
            cout << sdpdata.ele.bij[i][j] << " ";
        }
        cout << endl;
    }
    cout << " mat_info.coef" << endl;
    for(int j=0; j < sdpdata.ele.coef.size(); j++){
        cout << sdpdata.ele.coef[j] << " ";
    }
    cout << endl;
    cout << " mat_info.sup " << endl;
    sdpdata.ele.sup.disp();

    cout << "utnnz" << endl;
    for(int i=0; i < 3; i++){
        cout << "utnnz[" << i << "] = ";
        for(int j=0; j < sdpdata.utnnz[i].size(); j++){
            cout << sdpdata.utnnz[i][j] << " ";
        }
        cout << endl;
    }
}

void printMat(double *A, int ndim){
	cout << " A = " << endl;
	int pos = 0;
	for(int i=0; i<ndim; i++){
		for(int j=i; j<ndim; j++){
			printf("%4.2f ",A[pos]);
			pos = pos + 1;
		}
		printf("\n");
	}
	printf("\n");
}


void printOrder(vector<int> oriidx){
	cout << "ordering" << endl; 
	for(int i=0; i<oriidx.size(); i++){
		cout << oriidx[i] << " " ;
	}
	cout << endl;
}

void printExt(vector<vector<int> > extmat){
	printf("extmat\n");
	int ndim = extmat[0].size();
	for(int i=0; i<ndim; i++){
		for(int j=0; j<ndim; j++){
			printf("%d ", extmat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


void printCsp(vector<vector<double> > cspmat){
	printf("cspmat\n");
	int ndim = cspmat[0].size();
	for(int i=0; i<ndim; i++){
		for(int j=0; j<ndim; j++){
			printf("%7.4f ",cspmat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


void printScaleInfo(int slen, int blen, int mlen, vector<double> scalevalue, vector<double> bvect, vector<double> permmatrix){
	for(int i=0; i<slen; i++){
		cout << "scalevalue[" << i << "]= " << scalevalue[i] << endl;
	}
	for(int i=0; i<blen; i++){
		cout << "bvect[" << i << "]= " << bvect[i] << endl;
	}
	for(int i=0; i<mlen; i++){
		cout << "permmatrix[" << i << "]= " << permmatrix[i] << endl;
	}
}
