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
#include <cstdint>
#include "input.h"

void contraction(class polysystem & Polysys, vector<vector<double> > fixedVar, vector<int>& varMapping){
	int ndim = fixedVar[0].size();
	vector<int> idx;
	vector<double> val;
	for(int i=0; i<ndim; i++){
		//cout << fixedVar[0][i] << " ";
		if(fixedVar[0][i] == 0){
			idx.push_back(i);
		}
	}

	vector<poly>::iterator pit;
	vector<poly>::iterator pit_begin = Polysys.polynomial.begin();
	vector<poly>::iterator pit_end = Polysys.polynomial.end();
	list<mono>::iterator mit;
	list<mono>::iterator mit_begin;
	list<mono>::iterator mit_end;
	list<mono> monoList;
	class mono Mono;
	ndim = idx.size();
	int sizeCone;
	int pos = 0;
	for(pit = pit_begin; pit != pit_end; ++pit){
		sizeCone = (*pit).sizeCone;
		monoList = (*pit).monoList;
		mit_begin = monoList.begin();
		mit_end = monoList.end();
		(*pit).monoList.clear();
		(*pit).noTerms = 0;
		(*pit).dimVar = ndim;
		for(mit = mit_begin; mit != mit_end; ++mit){
			Mono.clear();
			Mono.allocCoef(sizeCone);
			Mono.allocSupp(ndim);
			for(int i=0; i<ndim; i++){
				//cout << "(nDim, idx[i], sup)  = (" << ndim << ", ";
				//cout << idx[i] << ", " << (*mit).Supp[idx[i]] << ")" << endl;
				Mono.setSupp(i,(*mit).Supp[idx[i]]);
			}
			for(int i=0; i<sizeCone; i++){
				Mono.setCoef((*mit).Coef[i],i);
			}
			(*pit).addMono(Mono);
		}
		//(*pit).setDegree();
		//(*pit).writePolyData();
	}
	Polysys.dimVar = ndim;
	vector<poly>::iterator pit2 = Polysys.polynomial.begin();
	int i=0;
	while(i < Polysys.numSys){ 
		if(i==0 && (*pit2).noTerms == 0){
			Mono.allocCoef(1);
			Mono.allocSupp(ndim);
			Mono.Coef[0] = 0;
			(*pit2).addMono(Mono);
			pos++;
			++pit2;
		}else if(i!=0 && (*pit2).noTerms == 0){
			Polysys.polynomial.erase(pit2);
			Polysys.numSys--;
			i--;	
		}else{
			//cout << pos << endl;
			Polysys.polynomial[i].noSys = pos;
			pos++;
			++pit2;
		}
		i++;
	}
	class bounds bds;
	bds.allocUpLo(ndim);
	for(int i=0; i<ndim; i++){
		//cout << "upper[" << i << "] = " << Polysys.bounds.ubd(idx[i]) << endl;
		//cout << "lower[" << i << "] = " << Polysys.bounds.lbd(idx[i]) << endl;
		bds.setUp(i+1,Polysys.bounds.ubd(idx[i]));	
		bds.setLow(i+1,Polysys.bounds.lbd(idx[i]));	
	}
	Polysys.bounds = bds;

	// Copy idx info to varMapping
	// Necessary so in conversion we know what variables have turned into what 
	varMapping = idx; 
}
void Substitute(class polysystem & Polysys, vector<int> idx,vector<double> val){
	int ndim = idx.size();
	vector<poly>::iterator pit;
	vector<poly>::iterator pit_begin = Polysys.polynomial.begin();
	vector<poly>::iterator pit_end = Polysys.polynomial.end();
	list<mono>::iterator mit;
	list<mono>::iterator mit_begin;
	list<mono>::iterator mit_end;
	list<mono> monoList;
	class poly Poly;	
	class mono Mono;	
	int power, sizeCone;
	
	//cout << "Substitute" << endl;

	for(int i=0; i<ndim; i++){
		for(pit = pit_begin; pit != pit_end; ++pit){
			sizeCone = (*pit).sizeCone;
			monoList = (*pit).monoList;
			mit_begin = monoList.begin();
			mit_end = monoList.end();
			(*pit).monoList.clear();
			(*pit).noTerms = 0;
			(*pit).degree = 0;
			for(mit = mit_begin; mit != mit_end; ++mit){
				power = (*mit).Supp[idx[i]];
				for(int k=0; k<sizeCone; k++){
					if(power > 0){
						(*mit).Coef[k] *= pow(val[i],power);
					}
				}
				(*mit).Supp[idx[i]] = 0;;
				(*pit).addMono((*mit));
			}
			(*pit).setDegree();
		} 
	}
}
void deleteVarMain(class polysystem & Polysys, vector<vector<double> > & fixedVar, bool & tf){
	// Step 1: find index of variables whose lower and upper bounds are same.
	//cout << "Step 1:" << endl;
	int ndim = fixedVar[0].size();
	vector<int> idx;
	vector<double> val;
	for(int i=0; i<ndim; i++){
		if(fixedVar[0][i] == 0 &&Polysys.bounds.ubd(i) - Polysys.bounds.lbd(i) < EPS){
			fixedVar[0][i] = 1;
			fixedVar[1][i] = Polysys.bounds.ubd(i);
			idx.push_back(i);
			val.push_back(Polysys.bounds.ubd(i));
		}
	}
	if(idx.empty() == false){
		Substitute(Polysys, idx, val);
		tf = true;
	}
	/*
	for(int i=0; i<idx.size(); i++){
		cout << "fixedVar : " << idx[i] << ", " << val[i] << endl;
	}
	*/
	// Step 2: find index of complementarity constraints.
	//cout << "Step 2:" << endl;
	idx.clear();
	val.clear();
	class sup Sup;
	double coef;
	int idx0, idx1;
	class bounds bds = Polysys.bounds;
	for(int i=1; i<Polysys.numSys;i++){
		//Polysys.polynomial[i].writePolyData();
		if(Polysys.polyIsComplimentarity(i,Sup,coef) == 1){
			if(fabs(coef) < EPS && Sup.nnz() == 2){
				idx0 = Sup.idx[0];
				idx1 = Sup.idx[1];
				if(fixedVar[0][idx1] == 0 && fixedVar[0][idx0] == 1 && fabs(fixedVar[1][idx0]) > EPS){
					fixedVar[0][idx1] = 1;
					fixedVar[1][idx1] = 0;
					idx.push_back(idx1);
					val.push_back(0);
				
				}
				if(fixedVar[0][idx0] == 0 && fixedVar[0][idx1] == 1 && fabs(fixedVar[1][idx1]) > EPS){
					fixedVar[0][idx0] = 1;
					fixedVar[1][idx0] = 0;
					idx.push_back(idx0);
					val.push_back(0);
				}


				if(fixedVar[0][idx0] == 0 && (bds.lbd(idx0) > 0 || bds.ubd(idx0) < 0)){
					fixedVar[0][idx0] = 1;
					fixedVar[1][idx0] = 0;
					idx.push_back(idx0);
					val.push_back(0);
				}
				if(fixedVar[0][idx1] == 0 && (bds.lbd(idx1) > 0 || bds.ubd(idx1) < 0)){
					fixedVar[0][idx1] = 1;
					fixedVar[1][idx1] = 0;
					idx.push_back(idx1);
					val.push_back(0);
				}
			}
		}
	}
	if(idx.empty() == false){
		Substitute(Polysys, idx, val);
		tf = true;
	}
	/*
	for(int i=0; i<idx.size(); i++){
		cout << "fixedVar : " << idx[i] << ", " << val[i] << endl;
	}
	*/
	// Step : find index of constraints x_i = constant.
	//cout << "Step 3:" << endl;
	idx.clear();
	val.clear();
	vector<int> comb;
	int noTerms, deg, typeCone;
	list<mono>::iterator mit0, mit1;
	list<mono> monoList;
	for(int i=1; i<Polysys.numSys;i++){
		//cout << i << endl;
		class mono Mono0, Mono1;
		class sup Sup0, Sup1;
		typeCone = Polysys.polynomial[i].typeCone;
		noTerms = Polysys.polynomial[i].noTerms;
		deg     = Polysys.polynomial[i].degree;
		monoList= Polysys.polynomial[i].monoList;
		//cout << "(noTerms, deg) = (" << noTerms << ", " << deg << ")" << endl;
		if(noTerms == 1 && deg == 1 && typeCone == -1){
			mit0 = monoList.begin();
			Mono0 = (*mit0);
			Mono0.getSuppComb(comb);
			//cout << "1: comb = " << comb[0] << endl;
			//Polysys.polynomial[i].writePolyData();
			if(fixedVar[0][comb[0]] ==0){
				fixedVar[0][comb[0]] = 1;
				fixedVar[1][comb[0]] = 0;
				idx.push_back(comb[0]);
				val.push_back(0);
			}
			comb.clear();	
		}else if(noTerms == 2 && deg == 1 && typeCone == -1){
			//Polysys.polynomial[i].writePolyData();
			mit0 = monoList.begin();
			Mono0 = *mit0;
			++mit0;
			mit1 = mit0;
			Mono1 = *mit1;
			Mono0.getSup(Sup0);
			Mono1.getSup(Sup1);
			//Sup0.disp();
			//Sup1.disp();
			//cout << "nnz0 = " << Sup0.nnz() << endl; 
			//cout << "nnz1 = " << Sup1.nnz() << endl; 
			if(Sup0.nnz() == 1 && Sup1.nnz() == 0){
				Mono0.getSuppComb(comb);
				//cout << "2: comb = " << comb[0] << endl;
				//Polysys.polynomial[i].writePolyData();
				//cout << comb[0] << endl;
				if(fixedVar[0][comb[0]] == 0){
					fixedVar[0][comb[0]] = 1;
					idx.push_back(comb[0]);
					if(fabs(Mono0.Coef[0]) > EPS){
						fixedVar[1][comb[0]] = -Mono1.Coef[0]/Mono0.Coef[0];
						val.push_back(-Mono1.Coef[0]/Mono0.Coef[0]);
					}else{
						cout << " ## Error: coef is sufficiently small."<< endl;
						exit(EXIT_FAILURE);
					}
				}
				comb.clear();
			}else if(Sup1.nnz() == 1 && Sup0.nnz() == 0){
				Mono1.getSuppComb(comb);
				//cout << "2: comb = " << comb[0] << endl;
				//Polysys.polynomial[i].writePolyData();
				//cout << comb[0] << endl;
				if(fixedVar[0][comb[0]] == 0){
					fixedVar[0][comb[0]] = 1;
					idx.push_back(comb[0]);
					if(fabs(Mono1.Coef[0]) > EPS){
						fixedVar[1][comb[0]] = -Mono0.Coef[0]/Mono1.Coef[0];
						val.push_back(-Mono0.Coef[0]/Mono1.Coef[0]);
					}else{
						cout << " ## Error: coef is sufficiently small."<< endl;
						exit(EXIT_FAILURE);
					}
				}
				comb.clear();
			}
		}
	}
	if(idx.empty() == false){
		Substitute(Polysys, idx, val);
		tf = true;
	}
/*
	for(int i=0; i<idx.size(); i++){
		cout << "fixedVar : " << idx[i] << ", " << val[i] << endl;
	}
*/	
}

void deleteVar(class polysystem & Polysys, vector<vector<double> > & fixedVar, vector<int>& varMapping){
	bool tf = true;
	while(tf == true){
		tf = false;
		deleteVarMain(Polysys, fixedVar, tf);
		//cout << " delete " << endl;
	}

	int nnz = 0;
	for(int i=0; i<fixedVar[0].size(); i++){
		if(fixedVar[0][i] == 1){
			nnz++;
		}
	}
	if(nnz > 0){
		cout << "# " << nnz << " variables fixed." << endl;
	}
		contraction(Polysys, fixedVar, varMapping);
}

void read_vec(std::ifstream& in, std::vector<int>& vec){
    std::uint64_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof n);
    vec.resize(static_cast<std::size_t>(n));
    in.read(reinterpret_cast<char*>(vec.data()), n * sizeof(int));
}

tuple<int,int, vector<int>, vector<int>, vector<int>> loadBlob(const string& filename){
    std::ifstream in(filename.c_str(), std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + filename);

    int a, b;
    std::vector<int> v1, v2, v3;

    in.read(reinterpret_cast<char*>(&a), sizeof a);
    in.read(reinterpret_cast<char*>(&b), sizeof b);

    read_vec(in, v1);
    read_vec(in, v2);
    read_vec(in, v3);
    return std::make_tuple(a, b, v1, v2, v3); 
}

void inputGMS(class polysystem & PolySys, string & fname){
	/* read gms part */

	string tmp,line,varLine,posLine,barLine,conLine,objFunc;
	int pos0,pos2,pos3,pos4,pos5,pos6,pos7,pos8,pos9,pos10;
	//int pos1;
	int pos11 = -1;
	int varnum,pvarnum,bvarnum,consnum, minOrmax;
	vector<string> var,pvar,bvar,cons,eq,lbd,ubd,fix;
	varnum = 0;
	pvarnum = 0;
	bvarnum = 0;
	consnum = 0;
	minOrmax = 1;
	string::iterator it;
	ifstream gms(fname.c_str());
	if(gms.is_open()){
		while(!gms.eof()){
			/* read line of input gms file */
			getline(gms,tmp);
			//pos1 = tmp.find("default levels");
			//if(pos1 > 0){
			//	break;
			//}
			/* find the end of line which SparsePOP can read */
			pos0 = tmp.find(";");
			/* find the comment part */	
			pos2 = tmp.find("*");	
			//cout << "pos2 = " << pos2 << endl;
			if( (pos2 >0|| pos2 == -1) && pos0 > 0){
				/* delete white spaces */
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),' '));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == ' '; });
				tmp.erase(it, tmp.end());	
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),'\n'));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == '\n'; });
				tmp = string(tmp.begin(), it);		
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),'\t'));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == '\t'; });
				tmp = string(tmp.begin(), it);		
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),'\r'));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == '\r'; });
				tmp = string(tmp.begin(), it);		
				if(tmp.size() > 0){
						line = line + tmp;
						//cout << "line = " << line << endl;
				}
				/* pick up the variables involved in POP */ 	
				pos3 = line.find("Variables");
				if(pos3 == -1){
					pos3 = line.find("variables");
				}
				pos4 = line.find("Positive");
				if(pos4 == -1){
					pos4 = line.find("positive");
				}
				pos5 = line.find("Binary");
				if(pos5 == -1){
					pos5 = line.find("binary");
				}
				/* read "Variables" */	
				if(pos3 >= 0 && pos4 == -1 && pos5 == -1){
					varLine = line.substr(pos3+9);
					read_vars(varLine,var,0);
					varnum = var.size();
				}
				/* read "Positive Variables" */
				if(pos3 >= 0 && pos4 >= 0 && pos5 == -1){
					posLine = line.substr(pos3+9);
					read_vars(posLine,pvar,1);
					pvarnum = pvar.size();
				}
				/* read "Binary Variables" */
				if(pos3 >= 0 && pos4 == -1 && pos5 >= 0){
					barLine = line.substr(pos3+9);
					read_vars(barLine,bvar,1);
					bvarnum = bvar.size();
				}
				/* read "Equations" */
				pos5 = line.find("Equations");
				if(pos5 >= 0){
					conLine = line.substr(pos5+9);
					read_vars(conLine,cons,1);
					consnum = cons.size();
				}
				/* find eq which has objective function. */
				/* read eq which has constraintsn. */
				pos6 = line.find("..");
				if(pos6 >= 0){
					conLine = line.substr(pos6+2);
					pos7 = conLine.find("objvar");
					if(pos7 >= 0){
						consnum = consnum-1;
						objFunc.assign(conLine);
					}else{
						eq.push_back(conLine);
					}
				}
				/* read lower and upper bounds */
				/* read fix bounds */
				pos8 = line.find("lo");
				pos9 = line.find("up");
				pos10 = line.find("fx");
				if(pos8 >= 0 && pos9 == -1 && pos10 == -1){
					lbd.push_back(line);
				}
				if(pos9 >= 0 && pos8 == -1 && pos10 == -1){
					ubd.push_back(line);
				}
				if(pos10 >= 0 && pos8 == -1 && pos9 == -1){
					fix.push_back(line);
				}
				if(pos8 != -1 && pos9 != -1){
					cout << " ## Error: bounds of input has something wrong." << endl;
					exit(EXIT_FAILURE);
				}else if(pos8 != -1 && pos10 != -1){
					cout << " ## Error: bounds of input has something wrong." << endl;
					exit(EXIT_FAILURE);

				}else if(pos9 != -1 && pos10 != -1){
					cout << " ## Error: bounds of input has something wrong." << endl;
					exit(EXIT_FAILURE);

				}
				if(pos11 == -1){
					pos11 = line.find("maximizingobjvar");
					if(pos11 != -1){
						minOrmax = -1;
					}else{
						minOrmax = 1;
					}
				}
				//cout << line << endl;
				//cout << pos11 << endl;
				//cout << "minOrmax = " << minOrmax << endl;
				line.clear();
			}else if((pos2 > 0|| pos2 == -1) && pos0 == -1){
				/* delete white spaces */
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),' '));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == ' '; });
				tmp.erase(it, tmp.end());	
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),'\n'));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == '\n'; });
				tmp.erase(it, tmp.end());	
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),'\t'));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == '\t'; });
				tmp.erase(it, tmp.end());	
				// it = remove_if(tmp.begin(), tmp.end(),bind2nd(equal_to<char>(),'\r'));
				it = remove_if(tmp.begin(), tmp.end(), [](char c){ return c == '\r'; });
				tmp.erase(it, tmp.end());
				if(tmp.size() > 0 ){
					//cout << " tmp.size = " << tmp.length() << endl;
					//cout << " tmp = " << tmp << endl;
					line = line + tmp;
					//cout << "line = " << line << endl;
				}	
				//	cout << "tmp = " << tmp << endl;
			}
		}
		if(varnum == 0){
			cout << " ## Error: Should input variables " << endl;
			exit(EXIT_FAILURE);
		}else if( consnum < eq.size()){
			cout << " ## Error: Some equations are not defined in Equations. " << endl;
			exit(EXIT_FAILURE);
		}else if(consnum > eq.size()){
			cout << " ## Error: Some equations are not defined in input file. "<< endl;
			exit(EXIT_FAILURE);
		}
	}else{
		cout << "Can not open file" << endl;
	}

	/* allocate POP from gms */
	int pos12, pos13, pos14;
	int scale = 1;	
	/* main */
	PolySys.numSys = 1+consnum+bvarnum;
	PolySys.dimVar = varnum;
	PolySys.polynomial.resize(1+consnum+bvarnum);
	/* allocate and set objective function */
	for(int i=0; i<1+consnum; i++){
		PolySys.polynomial[i].setNosysDimvar(i,varnum);
		vector<string> monos;
		if(i==0){
			detectParenthesis(objFunc);
			//cout << "objFunc = " << objFunc << endl;
			/* Separate objFunc into monomials */
			separate(objFunc, monos);
			PolySys.setMinOrMax(minOrmax);
			allocFromGMS(PolySys, i, scale,var,monos);
			/* Convert monomials into supports and coef*/
			/* allocate and set noTerms and degree from supports */
			PolySys.polynomial[i].setTypeSize(1,1);
		}else{
			scale = 0;
			detectParenthesis(eq[i-1]);
			separate(eq[i-1], monos);
			//allocFromGMS(PolySys, i, monos);
			pos14 = eq[i-1].find("=E=");
			pos12 = eq[i-1].find("=G=");
			pos13 = eq[i-1].find("=L=");
			if(pos14 >= 0){				
				PolySys.polynomial[i].setTypeSize(-1,1);
				scale = 1;
			}else if(pos12 >= 0){
				PolySys.polynomial[i].setTypeSize(1,1);
				scale = 1;
			}else if(pos13 >= 0){
				PolySys.polynomial[i].setTypeSize(1,1);
				scale = -1;
			}
			if(scale == 0){
				pos14 = eq[i-1].find("=e=");
				pos12 = eq[i-1].find("=g=");
				pos13 = eq[i-1].find("=l=");
				if(pos14 >= 0){
					PolySys.polynomial[i].setTypeSize(-1,1);
					scale = 1;
				}else if(pos12 >= 0){
					PolySys.polynomial[i].setTypeSize(1,1);
					scale = 1;
				}else if(pos13 >= 0){
					PolySys.polynomial[i].setTypeSize(1,1);
					scale = -1;
				}
			}
			if(scale == 0){
				cout << " ## Error: an equation does not have =E=, =L= and =G=." << endl;
				exit(EXIT_FAILURE);
			}
			allocFromGMS(PolySys, i, scale,var,monos);
		}
		PolySys.polynomial[i].setDegree();
	}
	//printStrVec(var);
	//printStrVec(bvar);
	set_binaryVars(PolySys, varnum, var, bvar);

	//PolySys.layawayObjConst();
	PolySys.posOflbds.resize(varnum);	
	PolySys.posOfubds.resize(varnum);	
	for(int i=0; i<varnum; i++){
		//PolySys.posOflbds[i].resize(varnum+1,0);
		//PolySys.posOfubds[i].resize(varnum+1,0);
		PolySys.posOflbds[i].resize(0);
		PolySys.posOfubds[i].resize(0);
	}	

	//set bounds
	set_bounds(PolySys,varnum,var,pvar,lbd,ubd,fix);
}

void allocMono(class mono & Mono, int scale, vector<string> var,string str){
	
	int size = var.size();
	Mono.allocSupp(size);
	Mono.allocCoef(1);
	int pow,length;		
	double value;
	vector<string> stvec;
	string tmp1,tmp2,tmp3;
	int pos, pos1, pos2, pos3;
	tmp1 = str;
	bool flag = false;
	while(true){
		pos = tmp1.find("*");
		length = tmp1.size();

		if(pos == -1){
			stvec.push_back(tmp1);
			break;
		}else if(pos == length-1){
			cout << str << endl;
			cout << " ## Error: input has something wrong." << endl;
			exit(EXIT_FAILURE);
		}else{
			tmp2 = tmp1.substr(0,pos);
			tmp1 = tmp1.substr(pos+1);
			stvec.push_back(tmp2);
		}
	}
	/*
	cout << "stvec --- " << endl;
	cout << "str = " << str << endl;
	printStrVec(stvec);	
	cout << "stvec --- " << endl;
	*/
	vector<string>::iterator it,varit;
	for(it = stvec.begin(); it != stvec.end(); ++it){
		pos1 = (*it).find("+");
		pos2 = (*it).find("-");
		if(pos1 == 0){
			tmp1 = (*it).substr(1);
			tmp2 = (*it).substr(0,1);
			/*
			cout << " tmp1= " << tmp1 << endl;
			cout << " tmp2= " << tmp2 << endl;
			*/
		}else if(pos2 == 0){
			tmp1 = (*it).substr(1);
			tmp2 = (*it).substr(0,1);
		}else{
			tmp1 = (*it);
			//cout << "tmp1 = " << tmp1 << endl;
		}
		pos3 = tmp1.find("^");
		if(pos3 != -1){
			//cout << "pos3 = " << pos3 << endl;
			tmp3 = tmp1.substr(pos3+1);
			tmp1 = tmp1.substr(0,pos3);
			//istringstream stream(tmp3);
			//stream >> pow;
			pow = atoi(tmp3.c_str());
			/*
			cout << *it << endl;
			cout << "tmp1 = " << tmp1 << endl;
			cout << "tmp3 = " << tmp3 << endl;
			cout << " pow = " << pow << endl;
			*/	
		}else{
			pow = 1;
		}
		pos3 = 0;
		/*
		cout << " *it = " << *it << endl;
		cout << " tmp1= " << tmp1 << endl;
		cout << " tmp2= " << tmp2 << endl;
		*/
		for(varit = var.begin(); varit != var.end(); ++varit){
			if(tmp1 == *varit){
				/*
				cout << "tmp1 = " << tmp1  << endl;
				cout << "varit= " << *varit<<endl;
				cout << "pow  = " << pow   <<endl;
				cout << "pos3 = " << pos3  <<endl;
				*/
				//Mono.writeMono();
				Mono.setSuppDense(pos3, pow);
				//Mono.Supp[pos3] = Mono.Supp[pos3] + pow;	
				//Mono.writeMono();
				flag = true;
				break;
			} 	
			pos3++;
		}
		double coef = Mono.Coef[0];
		if(fabs(coef) < EPS){
			coef = 1;
		}
		if(flag == false){
			if( (*it) == "+"){
				Mono.setCoef(+1*coef);
			}else if((*it) == "-"){
				Mono.setCoef(-1*coef);
			}else{
				//istringstream stream((*it));
				//stream >> value;
				value = atof((*it).c_str());
				Mono.setCoef(value*coef);
			}
		}else{
			if(tmp2 == "+"){
				Mono.setCoef(+1*coef);
			}else if(tmp2 == "-"){
				Mono.setCoef(-1*coef);
			}else{
				Mono.setCoef(+1*coef);
			}
		}
		flag = false;
		tmp1.clear();
		tmp2.clear();
		
	}
	Mono.convSuppDenseToSparse();
	/*	
	for(int i=0; i<Mono.supIdx.size(); i++){
		cout << Mono.supIdx[i] << " ";
		cout << Mono.supVal[i] << endl; 
	}
	Mono.writeMono();
	*/
	if(scale == -1){
		Mono.Coef[0] = scale*Mono.Coef[0];
	}	
	//Mono.writeMono();	
}

void allocFromGMS(class polysystem & PolySys, int i, int scale, vector<string> var, vector<string> monos){
	vector<string>::iterator it;
	class mono Mono;
	int pos1, pos2, pos3, val;
	if(i == 0){
		for(it = monos.begin(); it != monos.end(); ++it){
			pos1 = (*it).find("objvar");
			if(pos1 != -1){
				pos2 = (*it).find("-");
				if(pos2 != -1){
					scale = 1;
				}
				pos3 = (*it).find("+");
				if(pos2 == -1 && pos3 != -1){
					scale = -1;
				}
				monos.erase(it);
				break;
			}
		}
		val = PolySys.getMinOrMax();
		scale = val * scale;
	}

	for(it = monos.begin(); it != monos.end(); ++it){
		allocMono(Mono, scale, var, (*it));
		int numMonos = PolySys.polynomial[i].monoList.size();

		PolySys.polynomial[i].addMono(Mono);

		// For our purposes, this occurs when relations cancel due to flattening to the relational level
		// ex: Q(x,y) - Q(y,x) when flattened to the relational level becomes Q - Q which reduces to 0
		// So we must keep track of these cancellations so we can ground properly and use the information that is lost at the relational level
		if (PolySys.polynomial[i].monoList.size() < numMonos){
			// cout << "Yikes! Something cancelled in polynomial " << i << ", creating pair (" << Mono.supIdx[0] << "," << Mono.supVal[0] << ")" <<  endl;
			PolySys.polynomial[i].beenZero.push_back(make_pair(Mono.supIdx[0],Mono.supVal[0]));
			if (Mono.supIdx.size() > 1){
				cout << "Unhandled Case! Monomials of multiple atoms are cancelling in sparsePOP after reduction to the relational level" << endl;
			}
		}

		Mono.clear();
	}
}

void detectParenthesis(string str){
	int pos1,pos2;
	pos1 = str.find("(");
	pos2 = str.find(")");
	if(pos1 != -1 || pos2 != -1){
		cout << " ## Error: should delete parenthesis by your hand " << endl;
		exit(EXIT_FAILURE);	
	}
}


void separate(string str, vector<string> & monos){
	/* temporaly variables */
	string tmp,tmp1,tmp2,tmp3;
	int size,pos,pos1,pos2;
	double value;
	bool flag = true;	
	pos = str.find("=");
	if(pos == -1){
		cout << "line : " << str << endl;
		cout << "## Error: should set =E=, =G= or =L=." << endl;
		exit(EXIT_FAILURE);
	}else{
		size = str.length();
		tmp1 = str.substr(0,pos+3);
		//cout << "tmp1 = " << tmp1 << endl;
		tmp2 = str.substr(pos+3, size -(pos+3)-1);
		//cout << "tmp2 = " <<  tmp2 << endl;
		tmp3 = str.substr(0,pos);
		//cout << "tmp1=" << tmp1<< endl;
		//cout << "tmp2=" << tmp2<< endl;
		//cout << "tmp3=" << tmp3<< endl;
		if(tmp3 == "objvar"){
			tmp3 = str.substr(pos, 3);
			//cout << "tmp3=" << tmp3<< endl;
			tmp1 = tmp2 + tmp3;
			//cout << "tmp1=" << tmp1<< endl;
			//cout << "tmp2=" << tmp2<< endl;
		}else{
			value = atof(tmp2.c_str());
			if(fabs(value) > EPS){
				if(value > 0){
					tmp = tmp2.substr(0,1);
					if(tmp == "+"){
						tmp2 = tmp2.substr(1);
					}
					tmp2 = "-" + tmp2;
				}else{
					tmp2 = tmp2.substr(1);
					tmp2 = "+" + tmp2;
				}
				//cout << "tmp2 = " << tmp2 << endl; 
				monos.push_back(tmp2);
			}else{
				/*
				tmp = tmp2.substr(0,1);
				if(tmp ~= "-" && tmp ~= "+"){
					tmp2 = "+" + tmp2;
				}
				string plus = "+*";
				string minus= "-*";
				*/
			}
		}
	}	
	while(flag){
		pos1 = tmp1.find("+");
		pos2 = tmp1.find("-");
		if(pos1 > 0 && pos2 > 0){
			/* coef of the first monomial in eq does not have + */
			/* eq has + or - */
			pos = pos1;
			if(pos > pos2){
				pos = pos2;
			}
		}else if(pos1 > 0 && pos2 == -1){
			pos = tmp1.find("+",1);
		}else if(pos1 > 0 && pos2 == 0){
			pos2 = tmp1.find("-",1);
			pos = pos1;
			if(pos > pos2 && pos2 > 0){
				pos = pos2;
			}	
		}else if(pos1 == -1 && pos2 >0){
			/* coef of the first monomial in eq does not have + */
			/* eq has  - and does not have + */
			pos = pos2;
		}else if(pos1 == -1 && pos2 == -1){
			/* constraint consists of one monomial in the left hand side. */
			pos = tmp1.find("=");
		}else if(pos1 == 0){
			/* coef of the first monomial in eq has + */
			pos1 = tmp1.find("+",1);
			if(pos1 > 0){
				pos = pos1;
				if(pos > pos2 && pos2 >0){
					pos = pos2;
				}
			}else if(pos2 > 0){
				pos = pos2;
			}else{
				pos = tmp1.find("=");
			}
		}else if(pos2 == 0){
			/* coef of the first monomial in eq has - */
			/*
			cout << " tmp1 = " << tmp1 << endl;
			cout << " tmp2 = " << tmp2 << endl;
			*/
			pos2 = tmp1.find("-",1);
			if(pos2 > 0){
				pos = pos2;
				if(pos > pos1 && pos1 > 0){
					pos = pos1;
				}
			}else if(pos1 >0){
				pos = pos1;
			}else{
				pos = tmp1.find("=");
			}
		}
		tmp2 = tmp1.substr(0,pos);
		tmp1 = tmp1.substr(pos);
		monos.push_back(tmp2);
		/*
		cout << "tmp1 = " << tmp1 << endl;
		cout << "tmp2 = " << tmp2 << endl;
		*/
		pos = tmp1.find("=");
		if(pos == 0){
			flag = false;
		}	

	}
}


void printStrVec(vector<string> vec){
	vector<string>::iterator it;
	for(it= vec.begin(); it != vec.end(); ++it){
		cout << *it << endl;
	}
}
void read_vars(string line, vector<string> & vars, int i){
	int colon = line.find(";");
	line = line.substr(0,colon);
	//cout << line << endl;
	int comma = line.find(",");
	string tmp,newLine;
	newLine.assign(line);
	int size = line.length();
	//cout << " --- " << endl;
	while(comma >= 1){
		tmp.assign(newLine, 0,comma);
		newLine.assign(newLine,comma+1,size);
		size = newLine.length();
		if(tmp != "objvar"){
			vars.push_back(tmp);
		}
		//cout << "tmp = " << tmp << endl;
		comma = newLine.find(",");
	}
	/* pick up the last variable */
	/* need to drop colon */
	tmp.assign(newLine,0,size);
	//cout << "tmp = " << tmp << endl;
	if(i==0){
		if(tmp != "objvar"){
			vars.push_back(tmp);
			//cout << "tmp = " << tmp << endl;
		}
	}else{
		vars.push_back(tmp);
		//cout << "tmp = " << tmp << endl;
	}
	//cout << " --- " << endl;
	
}
void set_binaryVars(class polysystem & PolySys, int varnum, vector<string> var, vector<string> bvar){
	vector<string>::iterator ite;
	vector<int> idx;
	int bvarnum = bvar.size();
	idx.resize(bvarnum,0);
	int p = 0;
	//printStrVec(var);
	//printStrVec(bvar);
	//cout << "varnum = " << varnum << endl;
	for(int i=0; i<varnum; i++){
		for(ite = bvar.begin();ite !=bvar.end();++ite){
			if(var[i] == *ite){
				//cout << "(var, bvar) = (" << var[i]<<", " <<  *ite<< ")" << endl;
				idx[p] = i;
				p = p+1;
			}
		}
		if(p == bvarnum){
			break;
		}	
	}
	p = 0;
	//for(int i=0; i< idx.size(); i++){
	//	cout << idx[i] << endl;
	//}
	//cout << PolySys.numSys << endl;
	//cout << bvarnum << endl;
	for(int i=PolySys.numSys-bvarnum; i<PolySys.numSys; i++){
		PolySys.polynomial[i].setNoSys(i);
		PolySys.polynomial[i].setDimVar(varnum);
		PolySys.polynomial[i].setTypeSize(-1);
		PolySys.polynomial[i].setDegree(2);
		mono mono1, mono2;
		mono1.allocSupp(varnum);
		mono2.allocSupp(varnum);
		mono1.allocCoef(1);
		mono2.allocCoef(1);
		//cout << "varnum = " << varnum << endl;
		//cout << "idx[" << p << "] = " << idx[p] << endl;
		mono1.setSuppDense(idx[p], 1);
		mono2.setSuppDense(idx[p], 2);
		mono1.convSuppDenseToSparse();
		mono2.convSuppDenseToSparse();
		//mono1.setSupp(idx[p], 1);
		//mono2.setSupp(idx[p], 2);
		mono1.setCoef(-1);
		mono2.setCoef(1);
		//mono1.writeMono();
		//mono2.writeMono();
		PolySys.polynomial[i].addMono(mono1);
		PolySys.polynomial[i].addMono(mono2);
		mono1.clear();
		mono2.clear();
		p = p+1;
	}
	//PolySys.writePolynomials();
}

void set_bounds(class polysystem & PolySys, int varnum, vector<string> var, vector<string> pvar, vector<string> lbd, vector<string> ubd, vector<string> fix){

	//printStrVec(lbd);
	/* temporal variables */
	int pos;
	double value;
	string num, varName;
	vector<string>::iterator ite;

	PolySys.bounds.allocUpLo(varnum);
	PolySys.boundsNew.allocUpLo(varnum);
	for(int i=0; i<varnum; i++){
		//lower bounds
		for(ite = pvar.begin();ite !=pvar.end(); ++ite){
			if(var[i] == *ite){
				//cout << var[i] << "," << (*ite) << endl;
				PolySys.bounds.setLow(i+1,0);
				PolySys.boundsNew.setLow(i+1,0);
				break;
			}
		}
		for(ite = lbd.begin(); ite != lbd.end(); ++ite){
			pos = (*ite).find(".");
			/*
			cout << "length of str = " << (*ite).size() << endl;
			cout << " *ite   = " << *ite << endl;
			cout << " pos    = " << pos << endl;
			*/
			if(pos >= 0){
				varName = (*ite).substr(0,pos);
				//cout << "name = " << varName << endl;
				//cout << "var["<< i << "]  = " << var[i] << endl;
				if(var[i] == varName){
					pos     = (*ite).find("=");
					num = (*ite).substr(pos+1);
					pos = num.find(";");
					num = num.substr(0,pos);
					//cout << "Name = " << varName << endl;
					//cout << "num     = " << num << endl;
					//istringstream stream(num);
					//stream >> value;
					value = atof(num.c_str());
					if(value > PolySys.bounds.lbd(i)){
						PolySys.bounds.setLow(i+1,value);
						PolySys.boundsNew.setLow(i+1,value);
					}
					break;
				}
			}
		}
		//upper bounds
		for(ite = ubd.begin(); ite != ubd.end(); ++ite){
			pos = (*ite).find(".");
			if(pos >= 0){
				varName = (*ite).substr(0,pos);
				if(varName == var[i]){
					pos     = (*ite).find("=");
					num = (*ite).substr(pos+1);
					pos = num.find(";");
					num = num.substr(0,pos);
					//cout << "varName = " << varName << endl;
					//cout << "num     = " << num << endl;
					//istringstream stream(num);
					//stream >> value;
					value = atof(num.c_str());
					PolySys.bounds.setUp(i+1,value);
					PolySys.boundsNew.setUp(i+1,value);
					break;
				}
			}
		}
		//fixed case
		for(ite = fix.begin(); ite != fix.end(); ++ite){
			pos = (*ite).find(".");
			if(pos >= 0){
				varName = (*ite).substr(0,pos);
				if(varName == var[i]){
					pos     = (*ite).find("=");
					num = (*ite).substr(pos+1);
					pos = num.find(";");
					num = num.substr(0,pos);
					//cout << "varName = " << varName << endl;
					//cout << "num     = " << num << endl;
					//istringstream stream(num);
					//stream >> value;
					value = atoi(num.c_str());
					PolySys.bounds.setLow(i+1,value);
					PolySys.boundsNew.setLow(i+1,value);
					PolySys.bounds.setUp(i+1,value);
					PolySys.boundsNew.setUp(i+1,value);
					break;
				}
			}
		}
		// Output for bounds
		/*
		cout.precision(15);
		cout << "lbd[" << i << "]= " << PolySys.bounds.lbd(i) << endl;
		cout.precision(15);
		cout << "ubd[" << i << "]= " << PolySys.bounds.ubd(i) << endl;
		*/	
	}
}


/*	
	//parameter
	POP.param.relax_Order = 1;
	POP.param.multiCliquesFactor = 0;
	POP.param.eqTolerance = 0.0;
	POP.param.sparseSW = 0;
	POP.param.complementaritySW = 0;
	POP.param.boundSW = 2;
	POP.param.scalingSW = 1;
	POP.param.reduceMomentMatSW = 1;
	POP.param.perturbation = 0;
    POP.param.detailedInfFile = "";
*/

