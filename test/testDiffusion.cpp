#include "TestUtils.hpp"
#include "tug/Core/Matrix.hpp"
#include "gtest/gtest.h"
#include <gtest/gtest.h>
#include <stdexcept>
#include <tug/Diffusion.hpp>

#include <Eigen/src/Core/Matrix.h>
#include <string>

// include the configured header file
#include <testSimulation.hpp>

#define DIFFUSION_TEST(x) TEST(Diffusion, x)

using namespace Eigen;
using namespace std;
using namespace tug;

constexpr int row = 11;
constexpr int col = 11;

template <tug::APPROACH approach>
Diffusion<double, approach> setupSimulation(RowMajMat<double> &concentrations,
                                            double timestep, int iterations) {
  int domain_row = 10;
  int domain_col = 10;

  // Grid
  // RowMajMat<double> concentrations = MatrixXd::Constant(row, col, 0);
  concentrations(5, 5) = 1;

  Diffusion<double, approach> diffusiongrid(concentrations);

  diffusiongrid.getConcentrationMatrix() = concentrations;
  diffusiongrid.setDomain(domain_row, domain_col);

  diffusiongrid.setTimestep(timestep);
  diffusiongrid.setIterations(iterations);
  diffusiongrid.setDomain(domain_row, domain_col);

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
  diffusiongrid.setAlphaX(alpha);
  diffusiongrid.setAlphaY(alpha);

  return diffusiongrid;
}

constexpr double timestep = 0.001;
constexpr double iterations = 7000;

DIFFUSION_TEST(EqualityFTCS) {
  // set string from the header file
  string test_path = testSimulationCSVDir;
  RowMajMat<double> reference = CSV2Eigen(test_path);
  cout << "FTCS Test: " << endl;

  RowMajMat<double> concentrations = MatrixXd::Constant(row, col, 0);

  Diffusion<double, tug::FTCS_APPROACH> sim =
      setupSimulation<tug::FTCS_APPROACH>(concentrations, timestep, iterations);

  // Boundary bc = Boundary(grid);

  // Simulation

  // Diffusion<double, tug::FTCS_APPROACH> sim(grid, bc);
  // sim.setOutputConsole(CONSOLE_OUTPUT_ON);
  // sim.setTimestep(timestep);
  // sim.setIterations(iterations);
  sim.run();

  cout << endl;
  EXPECT_TRUE(checkSimilarity(reference, sim.getConcentrationMatrix(), 0.1));
}

DIFFUSION_TEST(EqualityBTCS) {
  // set string from the header file
  string test_path = testSimulationCSVDir;
  RowMajMat<double> reference = CSV2Eigen(test_path);
  cout << "BTCS Test: " << endl;

  RowMajMat<double> concentrations = MatrixXd::Constant(row, col, 0);

  Diffusion<double, tug::BTCS_APPROACH> sim =
      setupSimulation<tug::BTCS_APPROACH>(concentrations, timestep,
                                          iterations); // Boundary

  // Boundary bc = Boundary(grid);

  // Simulation
  // Diffusion<double, tug::FTCS_APPROACH> sim(grid, bc);
  // sim.setOutputConsole(CONSOLE_OUTPUT_ON);
  // sim.setTimestep(timestep);
  // sim.setIterations(iterations);
  sim.run();

  cout << endl;
  EXPECT_TRUE(checkSimilarityV2(reference, sim.getConcentrationMatrix(), 0.01));
}

DIFFUSION_TEST(InitializeEnvironment) {
  int rc = 12;
  RowMajMat<double> concentrations(rc, rc);
  // Grid64 grid(concentrations);
  // Boundary boundary(grid);

  EXPECT_NO_FATAL_FAILURE(Diffusion<double> sim(concentrations));
}

// DIFFUSION_TEST(SimulationEnvironment) {
//   int rc = 12;
//   Eigen::MatrixXd concentrations(rc, rc);
//   Grid64 grid(concentrations);
//   grid.initAlpha();
//   Boundary boundary(grid);
//   Diffusion<double, tug::FTCS_APPROACH> sim(grid, boundary);

//   EXPECT_EQ(sim.getIterations(), 1);

//   EXPECT_NO_THROW(sim.setIterations(2000));
//   EXPECT_EQ(sim.getIterations(), 2000);
//   EXPECT_THROW(sim.setIterations(-300), std::invalid_argument);

//   EXPECT_NO_THROW(sim.setTimestep(0.1));
//   EXPECT_DOUBLE_EQ(sim.getTimestep(), 0.1);
//   EXPECT_DEATH(sim.setTimestep(-0.3), ".* greater than zero.*");
// }

DIFFUSION_TEST(ClosedBoundaries) {
  constexpr std::uint32_t nrows = 5;
  constexpr std::uint32_t ncols = 5;

  RowMajMat<double> concentrations =
      RowMajMat<double>::Constant(nrows, ncols, 1.0);
  RowMajMat<double> alphax = RowMajMat<double>::Constant(nrows, ncols, 1E-5);
  RowMajMat<double> alphay = RowMajMat<double>::Constant(nrows, ncols, 1E-5);

  Diffusion<double> sim(concentrations);
  sim.getAlphaX() = alphax;
  sim.getAlphaY() = alphay;

  // tug::Grid64 grid(concentrations);

  // grid.setAlpha(alphax, alphay);

  // tug::Boundary bc(grid);
  auto &bc = sim.getBoundaryConditions();
  bc.setBoundarySideConstant(tug::BC_SIDE_LEFT, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_TOP, 1.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_BOTTOM, 1.0);

  // tug::Diffusion<double> sim(grid, bc);
  sim.setTimestep(1);
  sim.setIterations(1);

  RowMajMat<double> input_values(concentrations);
  sim.run();

  EXPECT_TRUE(
      checkSimilarityV2(input_values, sim.getConcentrationMatrix(), 1E-12));
}

DIFFUSION_TEST(ConstantInnerCell) {
  constexpr std::uint32_t nrows = 5;
  constexpr std::uint32_t ncols = 5;

  RowMajMat<double> concentrations =
      RowMajMat<double>::Constant(nrows, ncols, 1.0);
  RowMajMat<double> alphax = RowMajMat<double>::Constant(nrows, ncols, 1E-5);
  RowMajMat<double> alphay = RowMajMat<double>::Constant(nrows, ncols, 1E-5);

  Diffusion<double> sim(concentrations);
  sim.getAlphaX() = alphax;
  sim.getAlphaY() = alphay;

  // tug::Grid64 grid(concentrations);
  // grid.setAlpha(alphax, alphay);

  // tug::Boundary bc(grid);
  auto &bc = sim.getBoundaryConditions();
  // inner
  bc.setInnerBoundary(2, 2, 0);

  // tug::Diffusion<double> sim(grid, bc);
  sim.setTimestep(1);
  sim.setIterations(1);

  MatrixXd input_values(concentrations);
  sim.run();

  const auto &concentrations_result = sim.getConcentrationMatrix();

  EXPECT_DOUBLE_EQ(concentrations_result(2, 2), 0);
  EXPECT_LT(concentrations_result.sum(), input_values.sum());

  EXPECT_FALSE((concentrations_result.array() > 1.0).any());

  EXPECT_FALSE((concentrations_result.array() < 0.0).any());
}
