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
#include "global.h"
#include "conversion.h"
#include "polynomials.h"
#include "sup.h"
#include "spvec.h"
#include "input.h"
#include "cspGraph.h"

void printScaleInfo(int, int, int, vector<double>, vector<double>, vector<double>);
void printCsp(vector<vector<double> >);
void printExt(vector<vector<int> >);
void printOrder(vector<int>);
void printMat(double *, int);
void printSDPA(class mysdp);
void printDegOne(vector<int>);



