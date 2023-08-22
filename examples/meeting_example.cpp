#include "tug/Boundary.hpp"
#include <tug/Simulation.hpp>

int main(int argc, char *argv[]) {
    int row = 20;
    int col = 30;

    // grid
    Grid grid = Grid(row,col);
    MatrixXd concentrations = MatrixXd::Constant(row,col,0);
    concentrations(10,4) = 2000;
    grid.setConcentrations(concentrations);
    MatrixXd alphaX = MatrixXd::Constant(row,col,1);
    for (int i = 0; i < row; i++) {
        alphaX(i,0) = 0.05;
    }
    MatrixXd alphaY = MatrixXd::Constant(row,col,1);
    grid.setAlpha(alphaX, alphaY);

    // boundary
    Boundary bc = Boundary(grid);
    bc.setBoundarySideConstant(BC_SIDE_LEFT, 0);
    bc.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
    bc.setBoundarySideConstant(BC_SIDE_TOP, 0);
    bc.setBoundarySideConstant(BC_SIDE_BOTTOM, 0);
    bc.setBoundaryElementClosed(BC_SIDE_LEFT, 3);
    bc.setBoundaryElementClosed(BC_SIDE_LEFT, 4);
    bc.setBoundaryElementClosed(BC_SIDE_LEFT, 5);
    bc.setBoundaryElementConstant(BC_SIDE_BOTTOM, 17, 100);

    // simulation
    Simulation sim = Simulation(grid, bc, BTCS_APPROACH);
    sim.setTimestep(0.1);
    sim.setIterations(100);
    sim.setOutputCSV(CSV_OUTPUT_XTREME);

    sim.run();
}