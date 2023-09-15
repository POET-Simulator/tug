#include <Eigen/Eigen>
#include <tug/Simulation.hpp>

using namespace Eigen;

int main(int argc, char *argv[]) {
  int row = 20;
  int col = 20;
  Grid64 grid(row, col);

  MatrixXd concentrations = MatrixXd::Constant(row, col, 0);
  concentrations(10, 10) = 2000;
  grid.setConcentrations(concentrations);

  Boundary bc = Boundary(grid);
  bc.setBoundarySideClosed(BC_SIDE_LEFT);
  bc.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc.setBoundarySideClosed(BC_SIDE_TOP);
  bc.setBoundarySideClosed(BC_SIDE_BOTTOM);

  Simulation simulation = Simulation(grid, bc, CRANK_NICOLSON_APPROACH);
  simulation.setTimestep(0.1);
  simulation.setIterations(50);
  simulation.setOutputCSV(CSV_OUTPUT_XTREME);

  simulation.run();
}
