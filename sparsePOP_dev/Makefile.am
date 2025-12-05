bindir       = ${prefix}/bin
datadir      = ${prefix}/share/@PACKAGE_NAME@
exampledir   = ${prefix}/share/@PACKAGE_NAME@/example
docdir       = ${prefix}/share/doc/@PACKAGE_NAME@
data_DATA    = param.pop
example_DATA = example1.gms GMSfiles/*.gms
doc_DATA     = README PARAMETER.txt AUTHORS COPYING ChangeLog NEWS 

bin_PROGRAMS = sparsePOP
sparsePOP_SOURCES = clique.cpp conversion.cpp cspGraph.cpp debug.cpp info.cpp input.cpp mysdp.cpp polynomials.cpp sparsePOP.cpp spvec.cpp sup.cpp Parameters.cpp clique.h conversion.h cspGraph.h debug.h info.h  input.h mysdp.h polynomials.h spvec.h sup.h Parameters.h global.h

sparsePOP_LDADD = $(SDPA_LIB) $(MUMPS_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(METIS_LIB) $(COLAMD_LIB) $(LIBS)
sparsePOP_CXXFLAGS= $(METIS_INCLUDE) $(COLAMD_INCLUDE) $(SDPA_INCLUDE) $(MUMPS_INCLUDE) $(CXXFLAGS)



