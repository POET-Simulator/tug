#debian stable (currently bullseye)
cmake_minimum_required(VERSION 3.18)

project(tug CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)
find_package(easy_profiler)
option(EASY_OPTION_LOG "Verbose easy_profiler" 1)

## SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2 -mfma")
option(TUG_USE_OPENMP "Compile with OpenMP support" ON)

set(CMAKE_CXX_FLAGS_GENERICOPT "-O3 -march=native" CACHE STRING
 "Flags used by the C++ compiler during opt builds."
  FORCE)

set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
  "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel GenericOpt."
  FORCE)

option(TUG_USE_UNSAFE_MATH_OPT
  "Use compiler options to break IEEE compliances by
    oenabling reordering of instructions when adding/multiplying of floating
    points."
  OFF)

if(TUG_USE_UNSAFE_MATH_OPT)
  add_compile_options(-ffast-math)
endif()

option(TUG_ENABLE_TESTING
  "Run tests after succesfull compilation"
  ON)

add_subdirectory(src)

if(TUG_ENABLE_TESTING)
  add_subdirectory(test)
endif()

add_subdirectory(examples)
