cmake_minimum_required(VERSION 3.20)
project(sparsepop_test LANGUAGES CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE OFF)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# SparsePOP include paths
include_directories(
    ${CMAKE_SOURCE_DIR}  # Current directory (sparsePOP_dev)
    $ENV{HOME}/src/SuiteSparse/UFconfig
    $ENV{HOME}/src/SuiteSparse/AMD/Include
    $ENV{HOME}/src/SuiteSparse/COLAMD/Include
    $ENV{HOME}/src/MUMPS_4.10.0/include
    $ENV{HOME}/.local/include
)

# Library paths
link_directories(
    $ENV{HOME}/src/SuiteSparse/UFconfig
    $ENV{HOME}/src/SuiteSparse/AMD/Lib
    $ENV{HOME}/src/SuiteSparse/COLAMD/Lib
    $ENV{HOME}/src/MUMPS_4.10.0/lib
    $ENV{HOME}/.local/lib
)

# Build SparsePOP - note: no sparsePOP_dev/ prefix since we're already in that directory
add_executable(sparsePOP_test
    clique.cpp
    conversion.cpp
    cspGraph.cpp
    debug.cpp
    info.cpp
    input.cpp
    mysdp.cpp
    polynomials.cpp
    sparsePOP.cpp
    spvec.cpp
    sup.cpp
    Parameters.cpp
)

# Link libraries
target_link_libraries(sparsePOP_test
    sdpa           
    dmumps         
    mumps_common
    pord
    mpiseq
    metis      
    colamd        
    amd           
    openblas   
    lapack
    gfortran    
    pthread
    m
)

# Compiler flags
target_compile_options(sparsePOP_test PRIVATE
    -O2
    -funroll-all-loops
    -D_REENTRANT
)

# Runtime library path
set_target_properties(sparsePOP_test PROPERTIES
    INSTALL_RPATH "$ENV{HOME}/.local/lib:$ENV{HOME}/src/MUMPS_4.10.0/lib:$ENV{HOME}/src/SuiteSparse/COLAMD/Lib"
    BUILD_WITH_INSTALL_RPATH TRUE
)
