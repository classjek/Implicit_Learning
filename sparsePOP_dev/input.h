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
#ifndef __input_h__
#define __input_h__
#include "global.h"
#include "polynomials.h"
#include "conversion.h"
#include "Parameters.h"

void contraction(class polysystem &, vector<vector<double>>, vector<int>&);
void Substitute(class polysystem &, vector<int>, vector<double>);
void deleteVarMain(class polysystem &, vector<vector<double> > &, bool &);
void deleteVar(class polysystem &, vector<vector<double> > &, vector<int>&);
void read_vars(string, vector<string> &, int);
void printStrVec(vector<string>);
void set_bounds(class polysystem &, int, vector<string>, vector<string>, vector<string>, vector<string>, vector<string>);
void separate(string, vector<string> &);
void detectParenthesis(string);
void allocFromGMS(class polysystem &, int, int, vector<string>, vector<string>);
void allocMono(class mono &, int, vector<string>, string);

// for testing, can delete later
tuple<int,int, vector<int>, vector<int>, vector<int>> loadBlob(const string& filename = "forSparse.bin");

void inputGMS(class polysystem &,string &);
void readParam(class pop_params &, string, int);
void genCSP(class polysystem, vector<vector<double> > &);
void assignCmat(vector<vector<double> > &, vector<int>);
void intVecUnique(vector<int> &);
void set_binaryVars(class polysystem &, int, vector<string>, vector<string>);
#endif // __input_h__
