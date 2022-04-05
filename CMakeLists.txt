#debian stable (currently bullseye)
cmake_minimum_required(VERSION 3.18)

project(Diffusion CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)

option(USE_OPENMP "Compile with OpenMP support" ON)

set(CMAKE_CXX_FLAGS_GENERICOPT "-O3 -march=native -ffast-math" CACHE STRING
 "Flags used by the C++ compiler during opt builds."
  FORCE)

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel GenericOpt."
  FORCE)

add_subdirectory(app)
add_subdirectory(src)
