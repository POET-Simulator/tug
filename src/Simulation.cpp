#include <cmath>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <tug/Simulation.hpp>
#include <omp.h>

#include <fstream>

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "BTCSv2.cpp"


using namespace std;

Simulation::Simulation(Grid &grid, Boundary &bc, APPROACH approach) : grid(grid), bc(bc) {

    this->approach = approach;
    this->solver = THOMAS_ALGORITHM_SOLVER;
    this->timestep = -1; // error per default
    this->iterations = -1;
    this->innerIterations = 1;
    this->numThreads = omp_get_num_procs();
    
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

    if (approach == FTCS_APPROACH) {

        double deltaRowSquare;
        double deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();
        double minDeltaSquare;
        double maxAlphaX, maxAlphaY, maxAlpha;
        string dim;
        if (grid.getDim() == 2) {
            dim = "2D";

            deltaRowSquare = grid.getDeltaRow() * grid.getDeltaRow();

            minDeltaSquare = (deltaRowSquare < deltaColSquare) ? deltaRowSquare : deltaColSquare;
            maxAlphaX = grid.getAlphaX().maxCoeff();
            maxAlphaY = grid.getAlphaY().maxCoeff();
            maxAlpha = (maxAlphaX > maxAlphaY) ? maxAlphaX : maxAlphaY;

        } else if (grid.getDim() == 1) {
            dim = "1D";
            minDeltaSquare = deltaColSquare;
            maxAlpha = grid.getAlpha().maxCoeff();
            
        } else {
            throw_invalid_argument("Critical error: Undefined number of dimensions!");
        }

        // Courant-Friedrichs-Lewy condition
        double cfl = minDeltaSquare / (4*maxAlpha); 

        // stability equation from Wikipedia; might be useful if applied cfl does not work in some cases
        // double CFL_Wiki = 1 / (4 * maxAlpha * ((1/deltaRowSquare) + (1/deltaColSquare)));

        cout << "FTCS_" << dim << " :: CFL condition MDL: " << cfl << endl;
        cout << "FTCS_" << dim << " :: required dt=" << timestep <<  endl;

        if (timestep > cfl) {

            this->innerIterations = (int)ceil(timestep / cfl);
            this->timestep = timestep / (double)innerIterations;

            cerr << "Warning :: Timestep was adjusted, because of stability "
                    "conditions. Time duration was approximately preserved by "
                    "adjusting internal number of iterations."
                << endl;
            cout << "FTCS_" << dim << " :: Required " << this->innerIterations
                << " inner iterations with dt=" << this->timestep << endl;

        } else {

            this->timestep = timestep;
            cout << "FTCS_" << dim << " :: No inner iterations required, dt=" << timestep << endl;

        }

    } else {
        this->timestep = timestep;
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

void Simulation::setSolver(SOLVER solver) {
    if (this->approach == FTCS_APPROACH) {
        cerr << "Warning: Solver was set, but FTCS approach initialized. Setting the solver "
                "is thus without effect."
             << endl;
    }

    this->solver = solver;
}

void Simulation::setNumberThreads(int numThreads){
    if(numThreads>0 && numThreads<=omp_get_num_procs()){
        this->numThreads=numThreads;
    }
    else{
        int maxThreadNumber = omp_get_num_procs(); 
        
        throw_invalid_argument("Number of threads exceeds the number of processor cores or is less than 1.");
    }
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

    if (approach == FTCS_APPROACH) { // FTCS case

        for (int i = 0; i < iterations * innerIterations; i++) {
            if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                printConcentrationsConsole();
            }
            if (csv_output >= CSV_OUTPUT_VERBOSE) {
                printConcentrationsCSV(filename);
            }

            FTCS(this->grid, this->bc, this->timestep, this->numThreads);
    
            // if (i % (iterations * innerIterations / 100) == 0) {
            //     double percentage = (double)i / ((double)iterations * (double)innerIterations) * 100;
            //     if ((int)percentage % 10 == 0) {
            //         cout << "Progress: " << percentage << "%" << endl;
            //     }
            // }
        }

    } else if (approach == BTCS_APPROACH) { // BTCS case

        if (solver == EIGEN_LU_SOLVER) {
            for (int i = 0; i < iterations; i++) {
                if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                    printConcentrationsConsole();
                }
                if (csv_output >= CSV_OUTPUT_VERBOSE) {
                    printConcentrationsCSV(filename);
                }

                BTCS_LU(this->grid, this->bc, this->timestep, this->numThreads);
                
            }
        } else if (solver == THOMAS_ALGORITHM_SOLVER) {
            for (int i = 0; i < iterations; i++) {
                if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
                    printConcentrationsConsole();
                }
                if (csv_output >= CSV_OUTPUT_VERBOSE) {
                    printConcentrationsCSV(filename);
                }

                BTCS_Thomas(this->grid, this->bc, this->timestep, this->numThreads);
                
            }
        }

    } else if (approach == CRANK_NICOLSON_APPROACH) { // Crank-Nicolson case

        // TODO 

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

#endif
