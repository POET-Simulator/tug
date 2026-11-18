#include <Eigen/Eigen>
#include <tug/tug.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {
  // EASY_PROFILER_ENABLE;
  // profiler::startListen();
  // **************
  // **** GRID ****
  // **************
  // profiler::startListen();
  // create a grid with a 20 x 20 field
  int row = 40;
  int col = 50;

  // (optional) set the domain, e.g.:
  // grid.setDomain(20, 20);

  // (optional) set the concentrations, e.g.:
  // MatrixXd concentrations = MatrixXd::Constant(20,20,1000); //
  // #row,#col,value grid.setConcentrations(concentrations);
  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);
  concentrations(10, 10) = 2000;

  // ******************
  // **** BOUNDARY ****
  // ******************

  // set up a simulation environment
  Diffusion simulation(concentrations);

  // create a boundary with constant values
  auto &bc = simulation.getBoundaryConditions();
  bc.setBoundarySideClosed(BC_SIDE_LEFT);
  bc.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc.setBoundarySideClosed(BC_SIDE_TOP);
  bc.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // ************************
  // **** SIMULATION ENV ****
  // ************************

  // set the timestep of the simulation
  simulation.setTimestep(0.1); // timestep

  // set the number of iterations
  simulation.setIterations(300);

  // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT::XTREME);

  // **** RUN SIMULATION ****

  // run the simulation

  // EASY_BLOCK("SIMULATION")
  simulation.run();
  // EASY_END_BLOCK;
  // profiler::dumpBlocksToFile("test_profile.prof");
  // profiler::stopListen();
}
