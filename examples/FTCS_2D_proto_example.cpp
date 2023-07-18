#include <tug/Simulation.hpp>
#include <iostream>

int main(int argc, char *argv[]) {
    
    Grid grid = Grid(20,20);

    Boundary bc = Boundary(grid, BC_TYPE_CONSTANT);

    Simulation simulation = Simulation(grid, bc, FTCS_APPROACH);
    
    simulation.run();
    

}