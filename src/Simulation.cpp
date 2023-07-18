#include <stdexcept>
#include <tug/Simulation.hpp>

#include "FTCS.cpp"

using namespace std;

// auto FTCS(Grid &grid, Boundary &bc, double timestep) {

// }

Simulation::Simulation(Grid &grid, Boundary &bc, APPROACH approach) : grid(grid), bc(bc) {
    //probably to DEBUG assignment of grid and bc
    this->grid = grid;
    this->approach = approach;

    //TODO calculate max time step
    this->timestep = 0.01;
    this->iterations = 1000;
    this->csv_output = CSV_OUTPUT_OFF;
}

void Simulation::setOutputCSV(CSV_OUTPUT csv_output) {
    if (csv_output < CSV_OUTPUT_OFF && csv_output > CSV_OUTPUT_VERBOSE) {
        throw invalid_argument("Invalid CSV output option given!");
    }

    this->csv_output = csv_output;
}

void Simulation::setTimestep(double timestep) {
    //TODO check timestep in FTCS for max value
    this->timestep = timestep;
}

auto Simulation::getTimestep() {
    return this->timestep;
}

void Simulation::setIterations(int iterations) {
    this->iterations = iterations;
}

auto Simulation::getIterations() {
    return this->iterations;
}


auto Simulation::run() {
    if (approach == FTCS_APPROACH) {
        for (int i = 0; i < iterations; i++) {
            FTCS(grid, bc, timestep);
        }
    } else if (approach == BTCS_APPROACH) {
        for (int i = 0; i < iterations; i++) {
            //TODO
            break;
        }
    }
}
