#include <Eigen/Eigen>
#include <tug/tug.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {
  // **************
  // **** GRID ****
  // **************
  int row = 40;
  int col = 50;

  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);
  concentrations(10, 10) = 2000;

  // ************************
  // **** SIMULATION ENV ****
  // ************************
  Diffusion<double> simulation(concentrations);

  // ******************
  // **** BOUNDARY ****
  // ******************
  Boundary<double> &bc = simulation.getBoundaryConditions();
  bc.setBoundarySideClosed(BC_SIDE_LEFT);
  bc.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc.setBoundarySideClosed(BC_SIDE_TOP);
  bc.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // set the timestep of the simulation
  simulation.setTimestep(0.1);

  // set the number of iterations
  simulation.setIterations(300);

  // set kind of output [CSV_OUTPUT::OFF (default), CSV_OUTPUT::ON,
  // CSV_OUTPUT::VERBOSE, CSV_OUTPUT::XTREME]
  simulation.setOutputCSV(CSV_OUTPUT::XTREME);

  // **** RUN SIMULATION ****
  simulation.run();

  return 0;
}
