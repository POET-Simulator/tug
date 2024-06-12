/**
 * @file FTCS_2D_proto_example.cpp
 * @author Hannes Signer, Philipp Ungrund
 * @brief Creates a prototypical standard TUG simulation in 2D with FTCS
 * approach and constant boundary condition
 *
 */

#include <Eigen/Eigen>
#include <tug/Diffusion.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {

  // **************
  // **** GRID ****
  // **************

  // create a grid with a 20 x 20 field
  int row = 64;
  int col = 64;
  int n2 = row / 2 - 1;
  Grid64 grid(row, col);

  // (optional) set the domain, e.g.:
  // grid.setDomain(20, 20);

  // (optional) set the concentrations, e.g.:
  // MatrixXd concentrations = MatrixXd::Constant(20,20,1000); //
  // #row,#col,value
  MatrixXd concentrations = MatrixXd::Constant(row, col, 0);
  concentrations(n2, n2) = 1;
  concentrations(n2, n2 + 1) = 1;
  concentrations(n2 + 1, n2) = 1;
  concentrations(n2 + 1, n2 + 1) = 1;
  grid.setConcentrations(concentrations);

  // (optional) set alphas of the grid, e.g.:
  MatrixXd alphax = MatrixXd::Constant(row, col, 1E-4); // row,col,value
  MatrixXd alphay = MatrixXd::Constant(row, col, 1E-6); // row,col,value
  grid.setAlpha(alphax, alphay);

  // ******************
  // **** BOUNDARY ****
  // ******************

  // create a boundary with constant values
  Boundary bc = Boundary(grid);

  // (optional) set boundary condition values for one side, e.g.:
  bc.setBoundarySideConstant(BC_SIDE_LEFT, 0); // side,values
  bc.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
  bc.setBoundarySideConstant(BC_SIDE_TOP, 0);
  bc.setBoundarySideConstant(BC_SIDE_BOTTOM, 0);

  // ************************
  // **** SIMULATION ENV ****
  // ************************

  // set up a simulation environment
  Diffusion<double, tug::FTCS_APPROACH> simulation(
      grid, bc); // grid,boundary,simulation-approach

  // (optional) set the timestep of the simulation
  simulation.setTimestep(1000); // timestep

  // (optional) set the number of iterations
  simulation.setIterations(5);

  // (optional) set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT_OFF);

  // **** RUN SIMULATION ****

  // run the simulation
  simulation.run();

  return 0;
}
