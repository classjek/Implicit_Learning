C++ version of SparsePOP
(C++ Package of A Sparse SemiDefinite Programming Relaxation of 
Polynomial Optimization Problems)

Hayato Waki, Sunyoung Kim, Masakazu Kojima, 
Masakazu Muramatsu, Hiroshi Sugimoto and Makoto Yamashita

----------------------------------------
Index

1. Overview
2. System Requirements
3. Installation Guide
4. Usage of SparsePOP
5. Acknowledgements
6. License
7. E-mail address
8. History

----------------------------------------


1. Overview

This SparesPOP is a C++ package for a sparse semidefinite 
programming (SDP) relaxation method for approximating a 
global optimal value and solution of a polynomial optimization 
problem (POP) proposed by Waki et al. The sparse SDP relaxation 
exploits a sparse structure of polynomials in POPs when applying 
a hierarchy of LMI relaxations of increasing dimensions by 
Lasserre. The efficiency of SparsePOP to approximate optimal 
solutions of POPs is thus increased, and larger scale POPs can 
be handled. MATLAB is not required for this SparsePOP although
this is based on MATLAB version of SaprsePOP which has been 
already released.

1-1. Remark

For some POPs, the result by this SparsePOP may be different from 
one by MATLAB version even if you use symamd for both versions. 
The reason is that your COLAMD is different from COLAMD of 
MATLAB. See the follwoing site for more detailed:

http://www.cise.ufl.edu/research/sparse/colamd

Moreover, the following functions in MATLAB version of SparsePOP 
are not contained in this version, yet:
* functions for computing error bounds
* functions for refinement solutions by local optimization solvers.
* functions for checking linear independence of coefficinet matrix A.  


2. System Requirements

The following software packages are required for SparsePOP.

  i.   BLAS and LAPACK
	-- available from http://www.netlib.org/blas
	-- and from http://www.netlib.org/lapack
	-- we use dpptrf in liblapack.a for SparsePOP.	

  ii.  METIS 4.0, 4.0.1 or 4.0.3
	-- available from http://glaros.dtc.umn.edu/gkhome/metis/metis/overview 
	-- we use METIS_NodeND in libmetis.a for SparsePOP.
	-- One of references of METIS:
	``A Fast and Highly Quality Multilevel Scheme for Partitioning Irregular Graphs". 
	George Karypis and Vipin Kumar. 
	SIAM Journal on Scientific Computing, 
	Vol. 20, No. 1, pp. 359--392, 1999.

  iii. UFconfig
	-- available from http://www.cise.ufl.edu/research/sparse/UFconfig 
	-- we use UFconfig for COLAMD.

  iv. COLAMD 
	-- available from http://www.cise.ufl.edu/research/sparse/colamd
	-- we use symamd in libcolamd.a for SparsePOP. 

  v.  SDPA 7.3.1 or later
	-- available from http://sdpa.sourceforge.net/download.html
	-- we use sdpa for solving SDP relaxation problems.

NOTE1: If your OS is Debian or Ubuntu, we recommend to install UFconfig
and COLAMD by "apt-get". Also, because SuiteSparse contains them, you 
can install them by installing SuiteSparse. SuiteSparse is 
available from http://www.cise.ufl.edu/research/sparse/SuiteSparse/

NOTE2: SparsePOP can call some functions in METIS 4.0.3 or earlier. 
However, them in METIS 5 or later are not available in SparsePOP. 
We recommend to install METIS 4.0.3 or earlier. 


3. Installation Guide

Step 1: Execute ./configure 
To install SparsePOP, you need to give directories of libarary 
and include of all necesarry packages. If you do not give any 
optional package for ./configure, the following default values 
are set:

i.   --with-blas=/usr/lib 
     --with-lapack=/usr/lib
ii.  --with-metis-libdir=/usr/lib/ 
     --with-metis-includedir=/usr/include/metis
iii. --with-ufconfig=/usr/include/suitesparse
iv.  --with-colamd-libdir=/usr/lib 
     --with-colamd-includedir=/usr/include/suitesparse
v.   --with-sdpa-libdir=/usr/local/lib
     --with-sdpa-includedir=/usr/local/include
vi.  --with-mumps-sdpa-libdir=/usr/local/mumps/build/lib
     --with-mumps-sdpa-includedir=/usr/local/mumps/build/include

Step 2: Execute make
You can see an executable file sparsePOP in the directory.

3-1. Example of the installation

We assume that all necesarry packages are installed in the 
following directories:

i.   BLAS and LAPACK --> /usr/lib/
ii.  METIS --> /usr/local/lib/ 
               and /usr/local/include/
iii. UFconfig --> /usr/local/include/
iv.  COLAMD --> /usr/local/lib/
               and /usr/local/include/
v.   SDPA --> /usr/local/sdpa/lib/
              and /usr/local/sdpa/include/
vi.  MUMPS for SDPA --> /usr/local/sdpa/mumps/build/lib/
              and /usr/local/sdpa/mumps/build/include/

Then you can generate Makefile of SparsePOP by the following 
configure and its optional packages:

./configure --with-blas=/usr/lib --with-lapack=/usr/lib \
--with-metis-libdir=/usr/local/lib \
--with-metis-includedir=/usr/local/lib \
--with-ufconfig=/usr/local/include \
--with-colamd-libdir=/usr/local/lib \
--with-colamd-includedir=/usr/local/include \
--with-sdpa-libdir=/usr/local/sdpa/lib \
--with-sdpa-includedir=/usr/local/sdpa/incldue \
--with-mumps-sdpa-libdir=/usr/local/sdpa/mumps/build/lib \
--with-mumps-sdpa-includedir=/usr/local/sdpa/mumps/build/incldue

NOTE: If you use SDPA 7.3.1, then you do not need to add configure 
options for MUMPS of SDPA. Otherwise, user needs to add those 
configure options. 


4. Usage of SparsePOP

To execute SparsePOP, you need two parameter files. The first 
is `param.pop', which contains parameters for SDP relaxation 
and outputs of SparsePOP. The definitions of parameters in 
`param.pop' are given in `PARAMETER.txt'. The second is 
`param.sdpa', which is for sdpa. The detailed description is 
written in the SDPA manual. 

If you want to solve example1.gms, set relaxOrder = 3 and 
SDPsolverSW = 1 in `param.pop' and execute 

./sparsePOP example1.gms

Then SparsePOP displays infomation of SDP relaxation, the 
optimal value and a solution for example1.gms.


5. Acknowledgments

C++ version of SparsePOP was mainly developed 
while H. Waki stayed in CWI, Netherland. He thanks CWI 
and Professor Monique Laurent for the kind hospitality. 
The authors are grateful to Bartel Van Nieuwenhuyse
for reporting a bug in configure.in.

6. License

This software is distributed under GNU General Public License 2.0


7. E-mail address

kojima-spop@is.titech.ac.jp

Please send a message if you have any questions, ideas, or 
bug reports.


8. History

Version 2.59 (January 15, 2011)
* Add a reference of METIS in README.txt

* Some unsed variables in conversion.cpp, mysdp.cpp, etc. are eliminated.

Version 2.58 (December 15, 2011)
* Modified some functions in conversion.cpp, sup.cpp, 
polynomial.cpp and info.cpp. 

Version 2.56 (August 31, 2011)
* Fixed a small bug

Version 2.55 (July 30, 2011)
* SparsePOP can read the reserved keyword "maximizing objvar"
and "binary variables". If an input file contains the keyword 
"maximizing objvar", then SparsePOP deals with the POP as the 
maximization problem. In addition, if an input file contains 
the keyword "binary variables", then SparsePOP adds the polynomial 
equality x^2 -x = 0 in the original POP for such variables x.

* Some gams files in GLOBAL Library are added in the directory 
example/GMSformat. 

* For POPs which have constraints "x in {0,1}" and/or 
"x in {-1, +1}", SparsePOP can reduce the size of the 
resulting SDP relaxation problem by exploiting their 
constraints. This reduction is based on the following paper:

J.B.Lasserre, "An explicit equivalent positive semidefinite
program for nonlinear 0-1 programs", SIAM J.Optim., 12, 
pp. 756--769.

  For the reduction for a POP which has constraint "x in {0, 1}",
  set param.binarySW = 1. 
  For the reduction for a POP which has constraint "x in {-1, +1}",
  set param.SquareOneSW = 1. 
However, these reductions do not work when user computes error 
bounds of an approximated solution obtained by SparsePOP.

* We move "class pop_params" from conversion.h into Parameters.h.

Version 2.52 (June 21, 2011)
* Users can use multi threads for SDPA.
For example, if you want to use 2 threads, set OMP_NUM_THREADS as follows
before executing sparsePOP:
export OMP_NUM_THREADS=2
 
* Modify a part of outputting error messages.

Version 2.51 (June 2, 2011)
* Modify configure.in.
* Small bugs are fixed.

Version 2.50 (March 9, 2011)
* SparsePOPC++ released.
