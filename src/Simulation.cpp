#include <filesystem>
#include <stdexcept>
#include <string>
#include <tug/Simulation.hpp>

#include <fstream>

#include "FTCS.cpp"
#include "TugUtils.hpp"

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
    this->console_output = CONSOLE_OUTPUT_OFF;
    this->time_measure = TIME_MEASURE_OFF;
}

void Simulation::setOutputCSV(CSV_OUTPUT csv_output) {
    if (csv_output < CSV_OUTPUT_OFF && csv_output > CSV_OUTPUT_VERBOSE) {
        // throw invalid_argument("Invalid CSV output option given!");
        throw_invalid_argument("Invalid CSV output option given!");
    }

    this->csv_output = csv_output;
}

void Simulation::setOutputConsole(CONSOLE_OUTPUT console_output) {
    if (console_output < CONSOLE_OUTPUT_OFF && console_output > CONSOLE_OUTPUT_VERBOSE) {
        throw_invalid_argument("Invalid console output option given!");
    }

    this->console_output = console_output;
}

void Simulation::setTimeMeasure(TIME_MEASURE time_measure) {
    if (time_measure < TIME_MEASURE_OFF && time_measure > TIME_MEASURE_ON) {
        throw_invalid_argument("Invalid time measure option given!");
    }

    this->time_measure = time_measure;
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
    cout << grid.getConcentrations() << endl;
    cout << endl;
}

string Simulation::createCSVfile() {
    ofstream file;
    int appendIdent = 0;
    string appendIdentString;

     // APPROACH_ROW_COL_ITERATIONS
    string approachString = (approach == 0) ? "FTCS" : "BTCS";
    string row = to_string(grid.getRow());
    string col = to_string(grid.getCol());
    string numIterations = to_string(iterations);

    string filename = approachString + "_" + row + "_" + col + "_" + numIterations + ".csv";

    while (filesystem::exists(filename)) {
        appendIdent += 1;
        appendIdentString = to_string(appendIdent);
        filename = filename = approachString + "_" + row + "_" + col + "_" + numIterations + "-" + appendIdentString + ".csv";
    }

    file.open(filename);
    if (!file) {
        exit(1);
    }

    if (csv_output == CSV_OUTPUT_XTREME) {
        //rows
        //cols
        //iterations
        //boundary left
        //boundary right
        //boundary top
        //boundary bottom
        file << row << endl;
        file << col << endl;
        file << numIterations << endl;
        // TODO
        // file << to_string(bc.printBoundarySide) << endl;
        file << endl << endl;
    }

    file.close();

    return filename;
}

void Simulation::printConcentrationsCSV(string filename) {
    ofstream file;

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
    string filename;
    if (this->console_output > CONSOLE_OUTPUT_OFF) {
        printConcentrationsConsole();
    }
    if (this->csv_output > CSV_OUTPUT_OFF) {
        filename = createCSVfile();
        // printConcentrationsCSV(filename);
    }

    if (approach == FTCS_APPROACH) {
        
        for (int i = 0; i < iterations; i++) {
            if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsConsole();
            }
            if (csv_output == CSV_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsCSV(filename);
            }

            grid.setConcentrations(FTCS(grid, bc, timestep));
        }

    } else if (approach == BTCS_APPROACH) {

        for (int i = 0; i < iterations; i++) {
            if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsConsole();
            }
            if (csv_output == CSV_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsCSV(filename);
            }

            //TODO
            break;
        }

    }

    if (this->console_output > CONSOLE_OUTPUT_OFF) {
        printConcentrationsConsole();
    }
    if (this->csv_output > CSV_OUTPUT_OFF) {
        printConcentrationsCSV(filename);
    }
}
