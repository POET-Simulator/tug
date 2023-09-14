/**
 * @file FTCS_2D_proto_closed_mdl.cpp
 * @author Hannes Signer, Philipp Ungrund, MDL
 * @brief Creates a TUG simulation in 2D with FTCS approach and closed boundary
 * condition; optional command line argument: number of cols and rows
 *
 */

#include <cstdlib>
#include <iostream>
#include <tug/Simulation.hpp>

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

  // create a grid with a 20 x 20 field
  int n2 = row / 2 - 1;
  Grid grid = Grid(row, col);

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
  bc.setBoundarySideClosed(BC_SIDE_LEFT); // side,values
  bc.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc.setBoundarySideClosed(BC_SIDE_TOP);
  bc.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // ************************
  // **** SIMULATION ENV ****
  // ************************

  // set up a simulation environment
  Simulation simulation =
      Simulation(grid, bc, FTCS_APPROACH); // grid,boundary,simulation-approach

  // set the timestep of the simulation
  simulation.setTimestep(10000); // timestep

  // set the number of iterations
  simulation.setIterations(100);

  // (optional) set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT_VERBOSE);

  // **** RUN SIMULATION ****

  // run the simulation
  simulation.run();

  return 0;
}
