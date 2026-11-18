/**
 * @file FTCS_2D_proto_closed_mdl.cpp
 * @author Hannes Signer, Philipp Ungrund, MDL
 * @brief Creates a TUG simulation in 2D with FTCS approach and closed boundary
 * condition; optional command line argument: number of cols and rows
 *
 */

#include <Eigen/Eigen>
#include <cstdlib>
#include <iostream>
#include <tug/tug.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {

  int row = 64;

  if (argc == 2) {
    // no cmd line argument, take col=row=64
    row = atoi(argv[1]);
  }
  int col = row;

  std::cout << "Nrow =" << row << std::endl;
  // **************
  // **** GRID ****
  // **************
  int n2 = row / 2 - 1;

  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);
  concentrations(n2, n2) = 1;
  concentrations(n2, n2 + 1) = 1;
  concentrations(n2 + 1, n2) = 1;
  concentrations(n2 + 1, n2 + 1) = 1;

  // ************************
  // **** SIMULATION ENV ****
  // ************************
  Diffusion<double, FTCS_APPROACH> simulation(concentrations);

  RowMajMat<double> alphax = RowMajMat<double>::Constant(row, col, 1E-4);
  RowMajMat<double> alphay = RowMajMat<double>::Constant(row, col, 1E-6);
  simulation.setAlphaX(alphax);
  simulation.setAlphaY(alphay);

  // ******************
  // **** BOUNDARY ****
  // ******************
  Boundary<double> &bc = simulation.getBoundaryConditions();
  bc.setBoundarySideClosed(BC_SIDE_LEFT);
  bc.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc.setBoundarySideClosed(BC_SIDE_TOP);
  bc.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // set the timestep of the simulation
  simulation.setTimestep(10000);

  // set the number of iterations
  simulation.setIterations(100);

  // (optional) set kind of output [CSV_OUTPUT::OFF (default), CSV_OUTPUT::ON,
  // CSV_OUTPUT::VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT::VERBOSE);

  // **** RUN SIMULATION ****
  simulation.run();

  return 0;
}
