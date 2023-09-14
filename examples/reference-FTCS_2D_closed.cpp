#include "Eigen/Core"
#include <iostream>
#include <tug/Simulation.hpp>

using namespace std;

int main(int argc, char *argv[]) {
  int row = 50;
  int col = 50;
  int domain_row = 10;
  int domain_col = 10;

  // Grid
  Grid grid = Grid(row, col);
  grid.setDomain(domain_row, domain_col);

  MatrixXd concentrations = MatrixXd::Constant(row, col, 0);
  concentrations(5, 5) = 1;
  grid.setConcentrations(concentrations);

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
  grid.setAlpha(alpha, alpha);

  // Boundary
  Boundary bc = Boundary(grid);

  // Simulation
  Simulation sim = Simulation(grid, bc, FTCS_APPROACH);
  sim.setTimestep(0.001);
  sim.setIterations(10000);
  sim.setOutputCSV(CSV_OUTPUT_OFF);
  sim.setOutputConsole(CONSOLE_OUTPUT_OFF);

  // RUN
  sim.run();
}