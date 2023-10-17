/**
 * @file FTCS_2D_proto_example.cpp
 * @author Hannes Signer, Philipp Ungrund
 * @brief Creates a prototypical standard TUG simulation in 2D with FTCS
 * approach and constant boundary condition
 *
 */

#include <Eigen/Eigen>
#include <tug/Simulation.hpp>

using namespace Eigen;
using namespace tug;
// #include <easy/profiler.h>
// #define EASY_PROFILER_ENABLE ::profiler::setEnabled(true);

int main(int argc, char *argv[]) {
  // EASY_PROFILER_ENABLE;
  // profiler::startListen();
  // **************
  // **** GRID ****
  // **************
  // profiler::startListen();
  // create a grid with a 20 x 20 field
  int row = 20;
  int col = 20;
  Grid64 grid(row, col);

  // (optional) set the domain, e.g.:
  // grid.setDomain(20, 20);

  // (optional) set the concentrations, e.g.:
  // MatrixXd concentrations = MatrixXd::Constant(20,20,1000); //
  // #row,#col,value grid.setConcentrations(concentrations);
  MatrixXd concentrations = MatrixXd::Constant(row, col, 0);
  concentrations(0, 0) = 1999;
  grid.setConcentrations(concentrations);

  // (optional) set alphas of the grid, e.g.:
  // MatrixXd alphax = MatrixXd::Constant(20,20,1); // row,col,value
  // MatrixXd alphay = MatrixXd::Constant(20,20,1); // row,col,value
  // grid.setAlpha(alphax, alphay);

  // ******************
  // **** BOUNDARY ****
  // ******************

  // create a boundary with constant values
  Boundary bc = Boundary(grid);
  bc.setBoundarySideConstant(BC_SIDE_LEFT, 0);
  bc.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
  bc.setBoundarySideConstant(BC_SIDE_TOP, 0);
  bc.setBoundarySideConstant(BC_SIDE_BOTTOM, 0);

  // (optional) set boundary condition values for one side, e.g.:
  // VectorXd bc_left_values = VectorXd::Constant(20,1); // length,value
  // bc.setBoundaryConditionValue(BC_SIDE_LEFT, bc_left_values); // side,values
  // VectorXd bc_zero_values = VectorXd::Constant(20,0);
  // bc.setBoundaryConditionValue(BC_SIDE_LEFT, bc_zero_values);
  // bc.setBoundaryConditionValue(BC_SIDE_RIGHT, bc_zero_values);
  // VectorXd bc_front_values = VectorXd::Constant(20,2000);
  // bc.setBoundaryConditionValue(BC_SIDE_TOP, bc_front_values);
  // bc.setBoundaryConditionValue(BC_SIDE_BOTTOM, bc_zero_values);

  // ************************
  // **** SIMULATION ENV ****
  // ************************

  // set up a simulation environment
  Simulation simulation =
      Simulation<double, tug::FTCS_APPROACH>(grid, bc); // grid,boundary,simulation-approach

  // set the timestep of the simulation
  simulation.setTimestep(0.1); // timestep

  // set the number of iterations
  simulation.setIterations(10000);

  // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT_VERBOSE);

  // **** RUN SIMULATION ****

  // run the simulation

  // EASY_BLOCK("SIMULATION")
  simulation.run();
  // EASY_END_BLOCK;
  // profiler::dumpBlocksToFile("test_profile.prof");
  // profiler::stopListen();
}
