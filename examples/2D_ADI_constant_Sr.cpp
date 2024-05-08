#include "tug/Boundary.hpp"
#include "tug/Grid.hpp"
#include "tug/Simulation.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <cstdint>
#include <iostream>
int main() {

  constexpr std::uint32_t row = 100;
  constexpr std::uint32_t col = 200;

  tug::Grid64 grid(row, col);

  constexpr double sr_value = 4.5E-4;

  Eigen::MatrixXd concentrations =
      Eigen::MatrixXd::Constant(row, col, sr_value);
  grid.setConcentrations(concentrations);

  grid.setDomain(0.01, 0.02);

  constexpr double alpha_value = 1.1E-12;

  Eigen::MatrixXd alphax = Eigen::MatrixXd::Constant(row, col, alpha_value);
  Eigen::MatrixXd alphay = Eigen::MatrixXd::Constant(row, col, alpha_value);

  grid.setAlpha(alphax, alphay);

  tug::Boundary bc = tug::Boundary(grid);

  constexpr double sr_bc_value = 0.045;

  // north const
  for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(col / 2); i++) {
    bc.setBoundaryElementConstant(tug::BC_SIDE_TOP, i, sr_bc_value);
  }

  // south const
  for (std::uint32_t i = 1; i <= 19; i++) {
    bc.setBoundaryElementConstant(tug::BC_SIDE_BOTTOM, i, sr_bc_value);
  }

  // west const
  for (std::uint32_t i = 50; i <= 99; i++) {
    bc.setBoundaryElementConstant(tug::BC_SIDE_RIGHT, i, sr_bc_value);
  }

  tug::Simulation sim = tug::Simulation(grid, bc);
  sim.setTimestep(86400);
  sim.setIterations(20);

  sim.setOutputCSV(tug::CSV_OUTPUT_ON);

  sim.run();

  return 0;
}
