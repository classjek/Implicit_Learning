cmake_minimum_required(VERSION 3.20)
project(implicit_learning LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release or Debug" FORCE)
endif()

# OpenMP
find_package(OpenMP REQUIRED)

# Sources
set(SRC
  src/main.cpp
  src/executor.cpp
  src/domain.cpp
  src/kb_core.cpp
  src/observations.cpp
)

# Executable
add_executable(implicit_learning ${SRC})

# Include path now points to 'include' (not 'src')
target_include_directories(implicit_learning
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link OpenMP; also expose a compile-time define
target_link_libraries(implicit_learning PRIVATE OpenMP::OpenMP_CXX)
target_compile_definitions(implicit_learning PRIVATE HAVE_OPENMP=1)

# Warnings
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(implicit_learning PRIVATE -Wall -Wextra -Wpedantic)
elseif (MSVC)
  target_compile_options(implicit_learning PRIVATE /W4 /permissive-)
endif()


