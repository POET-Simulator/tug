#include "TestUtils.hpp"
#include <tug/Simulation.hpp>

#include <Eigen/src/Core/Matrix.h>
#include <doctest/doctest.h>
#include <stdio.h>
#include <string>

// include the configured header file
#include <testSimulation.hpp>

using namespace Eigen;
using namespace std;
using namespace tug;

Grid64 setupSimulation(double timestep, int iterations) {
  int row = 11;
  int col = 11;
  int domain_row = 10;
  int domain_col = 10;

  // Grid
  Grid grid = Grid64(row, col);
  grid.setDomain(domain_row, domain_col);

  MatrixXd concentrations = MatrixXd::Constant(row, col, 0);
  concentrations(5, 5) = 1;
  grid.setConcentrations(concentrations);

  MatrixXd alpha = MatrixXd::Constant(row, col, 1);
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 6; j++) {
      alpha(i, j) = 0.01;
    }
  }
  for (int i = 0; i < 5; i++) {
    for (int j = 6; j < 11; j++) {
      alpha(i, j) = 0.001;
    }
  }
  for (int i = 5; i < 11; i++) {
    for (int j = 6; j < 11; j++) {
      alpha(i, j) = 0.1;
    }
  }
  grid.setAlpha(alpha, alpha);

  return grid;
}

constexpr double timestep = 0.001;
constexpr double iterations = 7000;

TEST_CASE("equality to reference matrix with FTCS") {
  // set string from the header file
  string test_path = testSimulationCSVDir;
  MatrixXd reference = CSV2Eigen(test_path);
  cout << "FTCS Test: " << endl;

  Grid grid = setupSimulation(timestep, iterations); // Boundary
  Boundary bc = Boundary(grid);

  // Simulation

  Simulation sim = Simulation<double, tug::FTCS_APPROACH>(grid, bc);
  // sim.setOutputConsole(CONSOLE_OUTPUT_ON);
  sim.setTimestep(timestep);
  sim.setIterations(iterations);
  sim.run();

  cout << endl;
  CHECK(checkSimilarity(reference, grid.getConcentrations(), 0.1) == true);
}

TEST_CASE("equality to reference matrix with BTCS") {
  // set string from the header file
  string test_path = testSimulationCSVDir;
  MatrixXd reference = CSV2Eigen(test_path);
  cout << "BTCS Test: " << endl;

  Grid grid = setupSimulation(timestep, iterations); // Boundary
  Boundary bc = Boundary(grid);

  // Simulation
  Simulation sim = Simulation<double, tug::FTCS_APPROACH>(grid, bc);
  // sim.setOutputConsole(CONSOLE_OUTPUT_ON);
  sim.setTimestep(timestep);
  sim.setIterations(iterations);
  sim.run();

  cout << endl;
  CHECK(checkSimilarityV2(reference, grid.getConcentrations(), 0.01) == true);
}

TEST_CASE("Initialize environment") {
  int rc = 12;
  Grid64 grid(rc, rc);
  Boundary boundary(grid);

  CHECK_NOTHROW(Simulation sim(grid, boundary));
}

TEST_CASE("Simulation environment") {
  int rc = 12;
  Grid64 grid(rc, rc);
  Boundary boundary(grid);
  Simulation<double, tug::FTCS_APPROACH> sim(grid, boundary);

  SUBCASE("default paremeters") { CHECK_EQ(sim.getIterations(), -1); }

  SUBCASE("set iterations") {
    CHECK_NOTHROW(sim.setIterations(2000));
    CHECK_EQ(sim.getIterations(), 2000);
    CHECK_THROWS(sim.setIterations(-300));
  }

  SUBCASE("set timestep") {
    CHECK_NOTHROW(sim.setTimestep(0.1));
    CHECK_EQ(sim.getTimestep(), 0.1);
    CHECK_THROWS(sim.setTimestep(-0.3));
  }
}

TEST_CASE("Closed Boundaries - no change expected") {

  constexpr std::uint32_t nrows = 5;
  constexpr std::uint32_t ncols = 5;

  tug::Grid64 grid(nrows, ncols);

  auto concentrations = Eigen::MatrixXd::Constant(nrows, ncols, 1.0);
  auto alphax = Eigen::MatrixXd::Constant(nrows, ncols, 1E-5);
  auto alphay = Eigen::MatrixXd::Constant(nrows, ncols, 1E-5);

  grid.setConcentrations(concentrations);
  grid.setAlpha(alphax, alphay);

  tug::Boundary bc(grid);
  bc.setBoundarySideConstant(tug::BC_SIDE_LEFT, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_TOP, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_BOTTOM, 1.0);

  tug::Simulation<double> sim(grid, bc);
  sim.setTimestep(1);
  sim.setIterations(1);

  MatrixXd input_values(concentrations);
  sim.run();

  CHECK(checkSimilarityV2(input_values, grid.getConcentrations(), 1E-12) ==
        true);
}
