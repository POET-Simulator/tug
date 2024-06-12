find_package(doctest REQUIRED)

add_executable(testTug setup.cpp testDiffusion.cpp testGrid.cpp testFTCS.cpp testBoundary.cpp)
target_link_libraries(testTug doctest::doctest tug)

# get relative path of the CSV file
get_filename_component(testSimulationCSV "FTCS_11_11_7000.csv" REALPATH)
# set relative path in header file
configure_file(testSimulation.hpp.in testSimulation.hpp)
# include test directory with generated header file from above
target_include_directories(testTug PUBLIC "${CMAKE_CURRENT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/src")

add_custom_target(
  check
  COMMAND $<TARGET_FILE:testTug>
  DEPENDS testTug)
