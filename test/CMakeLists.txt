include(FetchContent)

FetchContent_Declare( 
  googletest 
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.15.2
)

FetchContent_MakeAvailable(googletest)


add_executable(testTug 
testSimulation.cpp 
testGrid.cpp 
testFTCS.cpp 
testBoundary.cpp
)
target_link_libraries(testTug tug GTest::gtest_main)
  
include(GoogleTest)
gtest_discover_tests(testTug)

# get relative path of the CSV file
get_filename_component(testSimulationCSV "FTCS_11_11_7000.csv" REALPATH)
# set relative path in header file
configure_file(testSimulation.hpp.in testSimulation.hpp)
# include test directory with generated header file from above
target_include_directories(testTug PUBLIC "${CMAKE_CURRENT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")