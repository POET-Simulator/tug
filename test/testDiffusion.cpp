#include "TestUtils.hpp"
#include <tug/Diffusion.hpp>

#include <Eigen/src/Core/Matrix.h>
#include <doctest/doctest.h>
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
  MatrixXd concentrations = MatrixXd::Constant(row, col, 0);
  concentrations(5, 5) = 1;

  Grid grid = Grid64(concentrations);
  grid.setDomain(domain_row, domain_col);

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

  Diffusion<double, tug::FTCS_APPROACH> sim(grid, bc);
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
  Diffusion<double, tug::FTCS_APPROACH> sim(grid, bc);
  // sim.setOutputConsole(CONSOLE_OUTPUT_ON);
  sim.setTimestep(timestep);
  sim.setIterations(iterations);
  sim.run();

  cout << endl;
  CHECK(checkSimilarityV2(reference, grid.getConcentrations(), 0.01) == true);
}

TEST_CASE("Initialize environment") {
  int rc = 12;
  Eigen::MatrixXd concentrations(rc, rc);
  Grid64 grid(concentrations);
  Boundary boundary(grid);

  CHECK_NOTHROW(Diffusion sim(grid, boundary));
}

TEST_CASE("Simulation environment") {
  int rc = 12;
  Eigen::MatrixXd concentrations(rc, rc);
  Grid64 grid(concentrations);
  grid.initAlpha();
  Boundary boundary(grid);
  Diffusion<double, tug::FTCS_APPROACH> sim(grid, boundary);

  SUBCASE("default paremeters") { CHECK_EQ(sim.getIterations(), 1); }

  SUBCASE("set iterations") {
    CHECK_NOTHROW(sim.setIterations(2000));
    CHECK_EQ(sim.getIterations(), 2000);
  }

  SUBCASE("set timestep") {
    CHECK_NOTHROW(sim.setTimestep(0.1));
    CHECK_EQ(sim.getTimestep(), 0.1);
  }
}

TEST_CASE("Closed Boundaries - no change expected") {

  constexpr std::uint32_t nrows = 5;
  constexpr std::uint32_t ncols = 5;

  auto concentrations = Eigen::MatrixXd::Constant(nrows, ncols, 1.0);
  auto alphax = Eigen::MatrixXd::Constant(nrows, ncols, 1E-5);
  auto alphay = Eigen::MatrixXd::Constant(nrows, ncols, 1E-5);

  tug::Grid64 grid(concentrations);

  grid.setAlpha(alphax, alphay);

  tug::Boundary bc(grid);
  bc.setBoundarySideConstant(tug::BC_SIDE_LEFT, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_TOP, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_BOTTOM, 1.0);

  tug::Diffusion<double> sim(grid, bc);
  sim.setTimestep(1);
  sim.setIterations(1);

  MatrixXd input_values(concentrations);
  sim.run();

  CHECK(checkSimilarityV2(input_values, grid.getConcentrations(), 1E-12) ==
        true);
}

TEST_CASE("Constant inner cell - 'absorbing' concentration") {
  constexpr std::uint32_t nrows = 5;
  constexpr std::uint32_t ncols = 5;

  auto concentrations = Eigen::MatrixXd::Constant(nrows, ncols, 1.0);
  auto alphax = Eigen::MatrixXd::Constant(nrows, ncols, 1E-5);
  auto alphay = Eigen::MatrixXd::Constant(nrows, ncols, 1E-5);

  tug::Grid64 grid(concentrations);
  grid.setAlpha(alphax, alphay);

  tug::Boundary bc(grid);
  // inner
  bc.setInnerBoundary(2, 2, 0);

  tug::Diffusion<double> sim(grid, bc);
  sim.setTimestep(1);
  sim.setIterations(1);

  MatrixXd input_values(concentrations);
  sim.run();

  CHECK(grid.getConcentrations()(2, 2) == 0);
  CHECK(grid.getConcentrations().sum() < input_values.sum());

  const bool greater_one = (grid.getConcentrations().array() > 1.0).any();
  CHECK(greater_one == false);

  const bool less_zero = (grid.getConcentrations().array() < 0.0).any();
  CHECK(less_zero == false);
}
