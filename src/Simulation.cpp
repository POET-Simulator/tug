#include <cmath>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <tug/Simulation.hpp>

#include <fstream>

#include "BTCSv2.cpp"
#include <tug/progressbar.hpp>


using namespace std;

Simulation::Simulation(Grid &grid, Boundary &bc, APPROACH approach) : grid(grid), bc(bc) {

    this->approach = approach;
    this->timestep = -1; // error per default
    this->iterations = -1;
    this->innerIterations = 1;
    
    this->csv_output = CSV_OUTPUT_OFF;
    this->console_output = CONSOLE_OUTPUT_OFF;
    this->time_measure = TIME_MEASURE_OFF;
}

void Simulation::setOutputCSV(CSV_OUTPUT csv_output) {
    if (csv_output < CSV_OUTPUT_OFF && csv_output > CSV_OUTPUT_VERBOSE) {
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
    if(timestep <= 0){
        throw_invalid_argument("Timestep has to be greater than zero.");
    }

    double deltaRowSquare;
    double deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();
    double minDeltaSquare;
    double maxAlphaX, maxAlphaY, maxAlpha;
    if (grid.getDim() == 2) {

        deltaRowSquare = grid.getDeltaRow() * grid.getDeltaRow();

        minDeltaSquare = (deltaRowSquare < deltaColSquare) ? deltaRowSquare : deltaColSquare;
        maxAlphaX = grid.getAlphaX().maxCoeff();
        maxAlphaY = grid.getAlphaY().maxCoeff();
        maxAlpha = (maxAlphaX > maxAlphaY) ? maxAlphaX : maxAlphaY;

    } else if (grid.getDim() == 1) {
        minDeltaSquare = deltaColSquare;
        maxAlpha = grid.getAlpha().maxCoeff();

        
    } else {
        throw_invalid_argument("Critical error: Undefined number of dimensions!");
    }

    // TODO check formula 1D case
    double CFL_MDL = minDeltaSquare / (4*maxAlpha); // Formula from Marco --> seems to be unstable
    double CFL_Wiki = 1 / (4 * maxAlpha * ((1/deltaRowSquare) + (1/deltaColSquare))); // Formula from Wikipedia

    cout << "FTCS_2D :: CFL condition MDL: " << CFL_MDL << endl;
    // cout << "FTCS_2D :: CFL condition Wiki: " << CFL_Wiki << endl;
    cout << "FTCS_2D :: required dt=" << timestep <<  endl;

    if (timestep > CFL_MDL) {

        this->innerIterations = (int)ceil(timestep / CFL_MDL);
        this->timestep = timestep / (double)innerIterations;

        cerr << "Warning: Timestep was adjusted, because of stability "
                "conditions. Time duration was approximately preserved by "
                "adjusting internal number of iterations."
             << endl;
        cout << "FTCS_2D :: Required " << this->innerIterations
            << " inner iterations with dt=" << this->timestep << endl;

    } else {

        this->timestep = timestep;
        cout << "FTCS_2D :: No inner iterations required, dt=" << timestep << endl;

    }

}

double Simulation::getTimestep() {
    return this->timestep;
}

void Simulation::setIterations(int iterations) {
    if(iterations <= 0){
        throw_invalid_argument("Number of iterations must be greater than zero.");
    }
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

    string approachString = (approach == 0) ? "FTCS" : "BTCS";
    string row = to_string(grid.getRow());
    string col = to_string(grid.getCol());
    string numIterations = to_string(iterations);

    string filename = approachString + "_" + row + "_" + col + "_" + numIterations + ".csv";

    while (filesystem::exists(filename)) {
        appendIdent += 1;
        appendIdentString = to_string(appendIdent);
        filename = approachString + "_" + row + "_" + col + "_" + numIterations + "-" + appendIdentString + ".csv";
    }

    file.open(filename);
    if (!file) {
        exit(1);
    }

    // adds lines at the beginning of verbose output csv that represent the boundary conditions and their values
    // -1 in case of closed boundary
    if (csv_output == CSV_OUTPUT_XTREME) {
        IOFormat one_row(StreamPrecision, DontAlignCols, "", " ");
        file << bc.getBoundarySideValues(BC_SIDE_LEFT).format(one_row) << endl; // boundary left
        file << bc.getBoundarySideValues(BC_SIDE_RIGHT).format(one_row) << endl; // boundary right
        file << bc.getBoundarySideValues(BC_SIDE_TOP).format(one_row) << endl; // boundary top
        file << bc.getBoundarySideValues(BC_SIDE_BOTTOM).format(one_row) << endl; // boundary bottom
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
    if (this->timestep == -1) {
        throw_invalid_argument("Timestep is not set!");
    }
    if (this->iterations == -1) {
        throw_invalid_argument("Number of iterations are not set!");
    }

    string filename;
    if (this->console_output > CONSOLE_OUTPUT_OFF) {
        printConcentrationsConsole();
    }
    if (this->csv_output > CSV_OUTPUT_OFF) {
        filename = createCSVfile();
    }

    auto begin = std::chrono::high_resolution_clock::now();

    if (approach == FTCS_APPROACH) {
        for (int i = 0; i < iterations * innerIterations; i++) {
            if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsConsole();
            }
            if (csv_output >= CSV_OUTPUT_VERBOSE) {
                printConcentrationsCSV(filename);
            }

            FTCS(this->grid, this->bc, this->timestep);
    
            if (i % (iterations * innerIterations / 100) == 0) {
                double percentage = (double)i / ((double)iterations * (double)innerIterations) * 100;
                if ((int)percentage % 10 == 0) {
                    cout << "Progress: " << percentage << "%" << endl;
                }
            }
        }

    } else if (approach == BTCS_APPROACH) {

        for (int i = 0; i < iterations * innerIterations; i++) {
            if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsConsole();
            }
            if (csv_output >= CSV_OUTPUT_VERBOSE) {
                printConcentrationsCSV(filename);
            }

            BTCS(this->grid, this->bc, this->timestep);
            
        }

    }

    auto end = std::chrono::high_resolution_clock::now();
    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    if (this->console_output > CONSOLE_OUTPUT_OFF) {
        printConcentrationsConsole();
    }
    if (this->csv_output > CSV_OUTPUT_OFF) {
        printConcentrationsCSV(filename);
    }
    if (this->time_measure > TIME_MEASURE_OFF) {
        string approachString = (approach == 0) ? "FTCS" : "BTCS";
        string dimString = (grid.getDim() == 1) ? "-1D" : "-2D";
        cout << approachString << dimString << ":: run() finished in " << milliseconds.count() << "ms" << endl;
    }
    
}
