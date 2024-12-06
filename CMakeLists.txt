# debian stable (currently bullseye)
cmake_minimum_required(VERSION 3.18)

project(
  tug
  VERSION 0.4
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(OpenMP)

include(GNUInstallDirs)

# find_package(easy_profiler) option(EASY_OPTION_LOG "Verbose easy_profiler" 1)

# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2 -mfma")
option(TUG_USE_OPENMP "Compile tug with OpenMP support" ON)

set(CMAKE_CXX_FLAGS_GENERICOPT
    "-O3 -march=native"
    CACHE STRING "Flags used by the C++ compiler during opt builds." FORCE)

set(CMAKE_BUILD_TYPE
    "${CMAKE_BUILD_TYPE}"
    CACHE
      STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel GenericOpt."
      FORCE)

option(
  TUG_USE_UNSAFE_MATH_OPT "Use compiler options to break IEEE compliances by
    oenabling reordering of instructions when adding/multiplying of floating
    points." OFF)


option(TUG_ENABLE_TESTING "Run tests after succesfull compilation" OFF)

option(TUG_HANNESPHILIPP_EXAMPLES "Compile example applications" OFF)

option(TUG_NAAICE_EXAMPLE "Enables NAAICE examples with export of CSV files"
       OFF)

add_library(tug INTERFACE)
target_include_directories(
  tug INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(tug INTERFACE Eigen3::Eigen)

target_compile_features(tug INTERFACE cxx_std_17)

if(TUG_USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(tug INTERFACE OpenMP::OpenMP_CXX)
endif()

if(TUG_USE_UNSAFE_MATH_OPT)
  target_compile_options(tug INTERFACE -ffast-math)
endif()

install(
  TARGETS tug
  EXPORT ${PROJECT_NAME}_Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "tugConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(
  EXPORT ${PROJECT_NAME}_Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/tug DESTINATION include)

if(TUG_ENABLE_TESTING)
enable_testing()
  add_subdirectory(test)
endif()

if(TUG_HANNESPHILIPP_EXAMPLES)
  add_subdirectory(examples)
endif()

if(TUG_NAAICE_EXAMPLE)
  add_subdirectory(naaice)
endif()
