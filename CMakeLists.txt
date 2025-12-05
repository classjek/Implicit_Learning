cmake_minimum_required(VERSION 3.20)
project(implicit_learning LANGUAGES CXX)

# C++ standard - CHANGED TO 17 for SparsePOP compatibility
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Position independent code OFF for SparsePOP
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release or Debug" FORCE)
endif()

# OpenMP
find_package(OpenMP REQUIRED)

# Your sources
set(MY_SRC
  src/main.cpp
  src/executor.cpp
  src/domain.cpp
  src/kb_core.cpp
  src/observations.cpp
  src/metrics.cpp
  src/spop.cpp
)

# SparsePOP sources - get all .cpp files from sparsePOP_dev/
file(GLOB SPARSEPOP_SRC "sparsePOP_dev/*.cpp")
# REMOVE SparsePOP's main.cpp (in SparsePOP.cpp, you commented it out, but this is extra safe)
list(REMOVE_ITEM SPARSEPOP_SRC "${CMAKE_CURRENT_SOURCE_DIR}/sparsePOP_dev/main.cpp")

# Executable with ALL sources (yours + SparsePOP)
add_executable(implicit_learning ${MY_SRC} ${SPARSEPOP_SRC})

# SparsePOP dependency paths
set(SDPA_DIR "/home/research/jakee/.local")
set(MUMPS_DIR "/home/research/jakee/src/MUMPS_4.10.0")
set(SUITESPARSE_DIR "/home/research/jakee/src/SuiteSparse")

# Include directories - ADD THESE FOR HEADERS
target_include_directories(implicit_learning PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/sparsePOP_dev
    ${SUITESPARSE_DIR}/UFconfig           # NEW
    ${SUITESPARSE_DIR}/AMD/Include        # NEW
    ${SUITESPARSE_DIR}/COLAMD/Include     # NEW
    ${MUMPS_DIR}/include                  # NEW
    ${SDPA_DIR}/include                   # NEW
)

# Link directories for SparsePOP dependencies
target_link_directories(implicit_learning PRIVATE
    ${SUITESPARSE_DIR}/UFconfig
    ${SUITESPARSE_DIR}/AMD/Lib
    ${SUITESPARSE_DIR}/COLAMD/Lib
    ${MUMPS_DIR}/lib
    ${SDPA_DIR}/lib
)

# Link all libraries (yours + SparsePOP dependencies)
target_link_libraries(implicit_learning PRIVATE
    OpenMP::OpenMP_CXX
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

# No PIE linker flag (required for MUMPS compatibility)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie")

# Compile definitions
target_compile_definitions(implicit_learning PRIVATE HAVE_OPENMP=1)

# Warnings
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(implicit_learning PRIVATE -Wall -Wextra -Wpedantic)
elseif (MSVC)
  target_compile_options(implicit_learning PRIVATE /W4 /permissive-)
endif()