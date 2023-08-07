#include <filesystem>
#include <stdexcept>
#include <string>
#include <tug/Simulation.hpp>

#include <fstream>

#include "FTCS.cpp"
#include "TugUtils.hpp"

using namespace std;

Simulation::Simulation(Grid &grid, Boundary &bc, APPROACH approach) : grid(grid), bc(bc) {

    this->approach = approach;


    // MDL no: we need to distinguish between "required dt" and
    // "number of (outer) iterations" at which the user needs an
    // output and the actual CFL-allowed timestep and consequently the
    // number of "inner" iterations which the explicit FTCS needs to
    // reach them. The following, at least at the moment, cannot be
    // computed here since "timestep" is not yet set when this
    // function is called. I brought everything into "FTCS_2D"!

    
    // TODO calculate max time step

    // double deltaRowSquare = grid.getDeltaRow() * grid.getDeltaRow();
    // double deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();

    // double minDelta2 = (deltaRowSquare < deltaColSquare) ? deltaRowSquare : deltaColSquare;
    // double maxAlphaX = grid.getAlphaX().maxCoeff();
    // double maxAlphaY = grid.getAlphaY().maxCoeff();
    // double maxAlpha = (maxAlphaX > maxAlphaY) ? maxAlphaX : maxAlphaY;

    // double CFL_MDL = minDelta2 / (4*maxAlpha); // Formula from Marco --> seems to be unstable
    // double CFL_Wiki = 1 / (4 * maxAlpha * ((1/deltaRowSquare) + (1/deltaColSquare))); // Formula from Wikipedia

    // cout << "Sim :: CFL condition MDL: " << CFL_MDL << endl;
    // double required_dt = this->timestep;
    // cout << "Sim :: required dt=" << required_dt <<  endl;
    // cout << "Sim :: CFL condition Wiki: " << CFL_Wiki << endl;

    // if (required_dt > CFL_MDL) {

    // 	int inner_iterations = (int) ceil(required_dt/CFL_MDL);
    // 	double allowed_dt = required_dt/(double)inner_iterations;
	    
    // 	cout << "Sim :: Required " << inner_iterations << " inner iterations with dt=" << allowed_dt <<  endl;
    // 	this->timestep = allowed_dt; 
    // 	this->iterations = inner_iterations;
    // } else {
    // 	cout << "Sim :: No inner iterations required, dt=" << required_dt <<  endl;
    // }
    
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
    if(timestep <= 0){
        throw_invalid_argument("Timestep has to be greater than zero.");
    }
    this->timestep = timestep;
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
    }

    if (approach == FTCS_APPROACH) {
        auto begin = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < iterations; i++) {
	    // MDL: distinguish between "outer" and "inner" iterations
	    std::cout << ":: run(): Outer iteration " << i+1 << "/" << iterations << endl;
            if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsConsole();
            }
            if (csv_output == CSV_OUTPUT_VERBOSE) {
                printConcentrationsCSV(filename);
            }

            FTCS(this->grid, this->bc, this->timestep);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

	// MDL: meaningful stdout messages
        std::cout << ":: run() finished in " << milliseconds.count() << "ms" << endl;

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
