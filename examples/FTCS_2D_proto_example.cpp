/**
 * @file FTCS_2D_proto_example.cpp
 * @author Hannes Signer, Philipp Ungrund
 * @brief Creates a prototypical standard TUG simulation in 2D with FTCS approach
 * and constant boundary condition
 * 
 */

#include <tug/Simulation.hpp>

int main(int argc, char *argv[]) {
    
    // **************
    // **** GRID ****
    // **************

    // create a grid with a 20 x 20 field
    int row = 20;
    int col = 20;
    Grid grid = Grid(row,col);

    // (optional) set the domain, e.g.:
    // grid.setDomain(20, 20);

    // (optional) set the concentrations, e.g.:
    // MatrixXd concentrations = MatrixXd::Constant(20,20,1000); // #row,#col,value
    // grid.setConcentrations(concentrations);

    // (optional) set alphas of the grid, e.g.:
    // MatrixXd alphax = MatrixXd::Constant(20,20,1); // row,col,value
    // MatrixXd alphay = MatrixXd::Constant(20,20,1); // row,col,value
    // grid.setAlpha(alphax, alphay);


    // ******************
    // **** BOUNDARY ****
    // ******************

    // create a boundary with constant values
    Boundary bc = Boundary(grid, BC_TYPE_CONSTANT);

    // (optional) set boundary condition values for one side, e.g.:
    // VectorXd bc_left_values = VectorXd::Constant(20,1); // length,value
    // bc.setBoundaryConditionValue(BC_SIDE_LEFT, bc_left_values); // side,values


    // ************************
    // **** SIMULATION ENV ****
    // ************************

    // set up a simulation environment
    Simulation simulation = Simulation(grid, bc, FTCS_APPROACH); // grid,boundary,simulation-approach

    // (optional) set the timestep of the simulation
    // simulation.setTimestep(0.01); // timestep

    // (optional) set the number of iterations
    simulation.setIterations(20);
    
    // **** RUN SIMULATION ****

    // run the simulation
    simulation.run();

}