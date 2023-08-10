#include <tug/Simulation.hpp>

int main(int argc, char *argv[]) {
    // **************
    // **** GRID ****
    // **************

    // create a linear grid with 20 cells
    int cells = 20;
    Grid grid = Grid(cells);

    MatrixXd concentrations = MatrixXd::Constant(1,20,0);
    concentrations(0,0) = 2000;
    // TODO add option to set concentrations with a vector in 1D case
    grid.setConcentrations(concentrations);


    // ******************
    // **** BOUNDARY ****
    // ******************

    // create a boundary with constant values
    Boundary bc = Boundary(grid);
    bc.setBoundarySideConstant(BC_SIDE_LEFT, 0);
    bc.setBoundarySideConstant(BC_SIDE_RIGHT, 0);


    // ************************
    // **** SIMULATION ENV ****
    // ************************

    // set up a simulation environment
    Simulation simulation = Simulation(grid, bc, BTCS_APPROACH); // grid,boundary,simulation-approach

    // set the timestep of the simulation
    simulation.setTimestep(0.1); // timestep

    // set the number of iterations
    simulation.setIterations(100);

    // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON, CSV_OUTPUT_VERBOSE]
    simulation.setOutputCSV(CSV_OUTPUT_VERBOSE);
    
    // **** RUN SIMULATION ****

    // run the simulation
    simulation.run();
}