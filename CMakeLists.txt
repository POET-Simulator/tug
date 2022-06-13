#debian stable (currently bullseye)
cmake_minimum_required(VERSION 3.18)

project(BTCSDiffusion CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)

## SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2 -mfma")
option(USE_OPENMP "Compile with OpenMP support" ON)

set(CMAKE_CXX_FLAGS_GENERICOPT "-O3 -march=native" CACHE STRING
 "Flags used by the C++ compiler during opt builds."
  FORCE)

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel GenericOpt."
  FORCE)

option(USE_UNSAFE_MATH_OPT
  "Use compiler options to break IEEE compliances by
    oenabling reordering of instructions when adding/multiplying of floating
    points."
  OFF)

if(USE_UNSAFE_MATH_OPT)
  add_compile_options(-ffast-math)
endif()

option(ENABLE_TESTING
  "Run tests after succesfull compilation"
  OFF)

add_subdirectory(app)
add_subdirectory(src)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
