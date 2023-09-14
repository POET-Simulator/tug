#include "TestUtils.cpp"
#include <doctest/doctest.h>
#include <stdio.h>
#include <string>
#include <tug/Simulation.hpp>

// include the configured header file
#include <testSimulation.hpp>

static Grid setupSimulation(APPROACH approach, double timestep,
                            int iterations) {
  int row = 11;
  int col = 11;
  int domain_row = 10;
  int domain_col = 10;

  // Grid
  Grid grid = Grid(row, col);
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

  // Boundary
  Boundary bc = Boundary(grid);

  // Simulation
  Simulation sim = Simulation(grid, bc, approach);
  // sim.setOutputConsole(CONSOLE_OUTPUT_ON);
  sim.setTimestep(timestep);
  sim.setIterations(iterations);
  sim.run();

  // RUN
  return grid;
}

TEST_CASE("equality to reference matrix with FTCS") {
  // set string from the header file
  string test_path = testSimulationCSVDir;
  MatrixXd reference = CSV2Eigen(test_path);
  cout << "FTCS Test: " << endl;
  Grid grid = setupSimulation(FTCS_APPROACH, 0.001, 7000);
  cout << endl;
  CHECK(checkSimilarity(reference, grid.getConcentrations(), 0.1) == true);
}

TEST_CASE("equality to reference matrix with BTCS") {
  // set string from the header file
  string test_path = testSimulationCSVDir;
  MatrixXd reference = CSV2Eigen(test_path);
  cout << "BTCS Test: " << endl;
  Grid grid = setupSimulation(BTCS_APPROACH, 1, 7);
  cout << endl;
  CHECK(checkSimilarityV2(reference, grid.getConcentrations(), 0.01) == true);
}

TEST_CASE("Initialize environment") {
  int rc = 12;
  Grid grid(rc, rc);
  Boundary boundary(grid);

  CHECK_NOTHROW(Simulation sim(grid, boundary, BTCS_APPROACH));
}

TEST_CASE("Simulation environment") {
  int rc = 12;
  Grid grid(rc, rc);
  Boundary boundary(grid);
  Simulation sim(grid, boundary, FTCS_APPROACH);

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
