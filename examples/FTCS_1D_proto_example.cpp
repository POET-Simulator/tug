#include <tug/Simulation.hpp>

int main(int argc, char *argv[]) {
    // **************
    // **** GRID ****
    // **************

    // create a linear grid with 20 cells
    int cells = 20;
    Grid grid = Grid(cells);

    MatrixXd concentrations = MatrixXd::Constant(1,20,20);
    grid.setConcentrations(concentrations);


    // ******************
    // **** BOUNDARY ****
    // ******************

    // create a boundary with constant values
    Boundary bc = Boundary(grid, BC_TYPE_CONSTANT);


    // ************************
    // **** SIMULATION ENV ****
    // ************************

    // set up a simulation environment
    Simulation simulation = Simulation(grid, bc, FTCS_APPROACH); // grid,boundary,simulation-approach

    // (optional) set the timestep of the simulation
    simulation.setTimestep(0.1); // timestep

    // (optional) set the number of iterations
    simulation.setIterations(100);

    // (optional) set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON, CSV_OUTPUT_VERBOSE]
    simulation.setOutputCSV(CSV_OUTPUT_OFF);
    
    // **** RUN SIMULATION ****

    // run the simulation
    simulation.run();
}