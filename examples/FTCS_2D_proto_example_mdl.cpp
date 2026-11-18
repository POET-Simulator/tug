/**
 * @file FTCS_2D_proto_example_mdl.cpp
 * @author Hannes Signer, Philipp Ungrund
 * @brief Creates a prototypical standard TUG simulation in 2D with FTCS
 * approach and constant boundary condition
 *
 */

#include <Eigen/Eigen>
#include <tug/tug.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {

  // **************
  // **** GRID ****
  // **************
  int row = 64;
  int col = 64;
  int n2 = row / 2 - 1;

  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);
  concentrations(n2, n2) = 1;
  concentrations(n2, n2 + 1) = 1;
  concentrations(n2 + 1, n2) = 1;
  concentrations(n2 + 1, n2 + 1) = 1;

  // ************************
  // **** SIMULATION ENV ****
  // ************************
  Diffusion<double, tug::FTCS_APPROACH> simulation(concentrations);

  RowMajMat<double> alphax = RowMajMat<double>::Constant(row, col, 1E-4);
  RowMajMat<double> alphay = RowMajMat<double>::Constant(row, col, 1E-6);
  simulation.setAlphaX(alphax);
  simulation.setAlphaY(alphay);

  // ******************
  // **** BOUNDARY ****
  // ******************
  Boundary<double> &bc = simulation.getBoundaryConditions();
  bc.setBoundarySideConstant(BC_SIDE_LEFT, 0);
  bc.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
  bc.setBoundarySideConstant(BC_SIDE_TOP, 0);
  bc.setBoundarySideConstant(BC_SIDE_BOTTOM, 0);

  // (optional) set the timestep of the simulation
  simulation.setTimestep(1000);

  // (optional) set the number of iterations
  simulation.setIterations(5);

  // (optional) set kind of output [CSV_OUTPUT::OFF (default), CSV_OUTPUT::ON,
  // CSV_OUTPUT::VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT::OFF);

  // **** RUN SIMULATION ****
  simulation.run();

  return 0;
}
