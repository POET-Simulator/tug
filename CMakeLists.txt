#debian stable (currently bullseye)
cmake_minimum_required(VERSION 3.18)

project(Diffusion CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)

## SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2 -mfma")
option(USE_OPENMP "Compile with OpenMP support" ON)

add_subdirectory(app)
add_subdirectory(src)
