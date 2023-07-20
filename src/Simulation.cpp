#include <stdexcept>
#include <tug/Simulation.hpp>

#include <fstream>

#include "FTCS.cpp"

using namespace std;

Simulation::Simulation(Grid grid, Boundary bc, APPROACH approach) : grid(grid), bc(bc) {
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

void Simulation::printConcentrationsConsole() {
    cout << "Concentrations:" << endl;
    cout << grid.getConcentrations() << endl;
    cout << endl;
}

void Simulation::printConcentrationsCSV(string ident) {
    ofstream file;

    string filename = "output-" + ident + ".csv";
    // string directory = "output/";
    file.open(filename, std::ios_base::app);
    if (!file) {
        exit(1);
    }

    IOFormat do_not_align(StreamPrecision, DontAlignCols);
    file << grid.getConcentrations().format(do_not_align) << endl;
    file << endl << endl;
    file.close();
}

void Simulation::run() {
    if (approach == FTCS_APPROACH) {
        printConcentrationsConsole();
        for (int i = 0; i < iterations; i++) {
            if (csv_output == CSV_OUTPUT_VERBOSE) {
                printConcentrationsCSV("test");
            }
            grid.setConcentrations(FTCS(grid, bc, timestep));
            // if (i != 0 && i % 200 == 0) {
            //     printConcentrationsConsole();
            // }
        }
        printConcentrationsConsole();
        if (csv_output >= CSV_OUTPUT_ON) {
            printConcentrationsCSV("test");
        }
    } else if (approach == BTCS_APPROACH) {
        for (int i = 0; i < iterations; i++) {
            //TODO
            break;
        }
    }
}
