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

    double deltaRowSquare = grid.getDeltaRow() * grid.getDeltaRow();
    double deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();

    double minDelta = (deltaRowSquare < deltaColSquare) ? deltaRowSquare : deltaColSquare;
    double maxAlphaX = grid.getAlphaX().maxCoeff();
    double maxAlphaY = grid.getAlphaY().maxCoeff();
    double maxAlpha = (maxAlphaX > maxAlphaY) ? maxAlphaX : maxAlphaY;

    //double maxStableTimestep = minDelta / (2*maxAlpha); // Formula from Marco --> seems to be unstable
    double maxStableTimestep = 1 / (4 * maxAlpha * ((1/deltaRowSquare) + (1/deltaColSquare))); // Formula from Wikipedia

    cout << maxStableTimestep << endl;

    this->timestep = maxStableTimestep;

    
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

double Simulation::getTimestep() {
    return this->timestep;
}

void Simulation::setIterations(int iterations) {
    this->iterations = iterations;
}

int Simulation::getIterations() {
    return this->iterations;
}

void Simulation::printConcentrationsConsole() {
    cout << "Concentrations:" << endl;
    cout << grid.getConcentrations() << endl;
    cout << endl;
}

void Simulation::printConcentrationsCSV(string ident, bool appendMode) {
    ofstream file;

    string filename = "output-" + ident + ".csv";
    // string directory = "output/";

    if (appendMode) {
        file.open(filename, std::ios_base::app);
    } else {
        file.open(filename);
    }
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
                printConcentrationsCSV("test", true);
            }
            grid.setConcentrations(FTCS(grid, bc, timestep));
            // if (i != 0 && i % 200 == 0) {
            //     printConcentrationsConsole();
            // }
        }
        printConcentrationsConsole();
        if (csv_output >= CSV_OUTPUT_ON) {
            bool append = false;
            if (csv_output == CSV_OUTPUT_VERBOSE) {
                append = true;
            }
            printConcentrationsCSV("test", append);
        }
    } else if (approach == BTCS_APPROACH) {
        for (int i = 0; i < iterations; i++) {
            //TODO
            break;
        }
    }
}
