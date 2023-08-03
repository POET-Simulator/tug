find_library(DOCTEST_LIB doctest)

if(NOT DOCTEST_LIB)
  include(FetchContent)

  FetchContent_Declare(
    DocTest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.9)

  FetchContent_MakeAvailable(DocTest)
endif()

add_executable(testTug setup.cpp testSimulation.cpp) # testBoundaryCondition.cpp testDiffusion.cpp
target_link_libraries(testTug doctest tug)

# get relative path of the CSV file
get_filename_component(testSimulationCSV "FTCS_11_11_7000.csv" REALPATH CACHE)
# set relative path in header file
configure_file(testSimulation.hpp.in testSimulation.hpp)
# include test directory with generated header file from above
target_include_directories(testTug PUBLIC "${CMAKE_CURRENT_BINARY_DIR}")

add_custom_target(
  check
  COMMAND $<TARGET_FILE:testTug>
  DEPENDS testTug)
