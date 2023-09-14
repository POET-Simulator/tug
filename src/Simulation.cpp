#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <tug/Simulation.hpp>

#include "Schemes.hpp"
#include "TugUtils.hpp"

Simulation::Simulation(Grid &_grid, Boundary &_bc, APPROACH _approach)
    : grid(_grid), bc(_bc), approach(_approach) {}

void Simulation::setOutputCSV(CSV_OUTPUT csv_output) {
  if (csv_output < CSV_OUTPUT_OFF && csv_output > CSV_OUTPUT_VERBOSE) {
    throw_invalid_argument("Invalid CSV output option given!");
  }

  this->csv_output = csv_output;
}

void Simulation::setOutputConsole(CONSOLE_OUTPUT console_output) {
  if (console_output < CONSOLE_OUTPUT_OFF &&
      console_output > CONSOLE_OUTPUT_VERBOSE) {
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
  if (timestep <= 0) {
    throw_invalid_argument("Timestep has to be greater than zero.");
  }

  if (approach == FTCS_APPROACH || approach == CRANK_NICOLSON_APPROACH) {

    const double deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();
    // will be 0 if 1D, else ...
    const double deltaRowSquare = grid.getDim() != 1
                                      ? grid.getDeltaRow() * grid.getDeltaRow()
                                      : deltaColSquare;
    const double minDeltaSquare =
        (deltaRowSquare < deltaColSquare) ? deltaRowSquare : deltaColSquare;

    double maxAlpha = std::numeric_limits<double>::quiet_NaN();

    // determine maximum alpha
    if (grid.getDim() == 2) {

      const double maxAlphaX = grid.getAlphaX().maxCoeff();
      const double maxAlphaY = grid.getAlphaY().maxCoeff();
      maxAlpha = (maxAlphaX > maxAlphaY) ? maxAlphaX : maxAlphaY;

    } else if (grid.getDim() == 1) {
      maxAlpha = grid.getAlpha().maxCoeff();

    } else {
      throw_invalid_argument("Critical error: Undefined number of dimensions!");
    }

    const std::string dim = std::to_string(grid.getDim()) + "D";

    // Courant-Friedrichs-Lewy condition
    double cfl = minDeltaSquare / (4 * maxAlpha);

    // stability equation from Wikipedia; might be useful if applied cfl does
    // not work in some cases double CFL_Wiki = 1 / (4 * maxAlpha *
    // ((1/deltaRowSquare) + (1/deltaColSquare)));

    const std::string &approachPrefix = this->approach_names[approach];
    std::cout << approachPrefix << "_" << dim << " :: CFL condition: " << cfl
              << std::endl;
    std::cout << approachPrefix << "_" << dim << " :: required dt=" << timestep
              << std::endl;

    if (timestep > cfl) {

      this->innerIterations = (int)ceil(timestep / cfl);
      this->timestep = timestep / (double)innerIterations;

      std::cerr << "Warning :: Timestep was adjusted, because of stability "
                   "conditions. Time duration was approximately preserved by "
                   "adjusting internal number of iterations."
                << std::endl;
      std::cout << approachPrefix << "_" << dim << " :: Required "
                << this->innerIterations
                << " inner iterations with dt=" << this->timestep << std::endl;

    } else {

      this->timestep = timestep;
      std::cout << approachPrefix << "_" << dim
                << " :: No inner iterations required, dt=" << timestep
                << std::endl;
    }

  } else {
    this->timestep = timestep;
  }
}

double Simulation::getTimestep() { return this->timestep; }

void Simulation::setIterations(int iterations) {
  if (iterations <= 0) {
    throw_invalid_argument("Number of iterations must be greater than zero.");
  }
  this->iterations = iterations;
}

void Simulation::setSolver(SOLVER solver) {
  if (this->approach == FTCS_APPROACH) {
    std::cerr
        << "Warning: Solver was set, but FTCS approach initialized. Setting "
           "the solver "
           "is thus without effect."
        << std::endl;
  }

  this->solver = solver;
}

void Simulation::setNumberThreads(int numThreads) {
  if (numThreads > 0 && numThreads <= omp_get_num_procs()) {
    this->numThreads = numThreads;
  } else {
    int maxThreadNumber = omp_get_num_procs();
    std::string outputMessage =
        "Number of threads exceeds the number of processor cores (" +
        std::to_string(maxThreadNumber) + ") or is less than 1.";

    throw_invalid_argument(outputMessage);
  }
}

int Simulation::getIterations() { return this->iterations; }

void Simulation::printConcentrationsConsole() {
  std::cout << grid.getConcentrations() << std::endl;
  std::cout << std::endl;
}

std::string Simulation::createCSVfile() {
  std::ofstream file;
  int appendIdent = 0;
  std::string appendIdentString;

  // string approachString = (approach == 0) ? "FTCS" : "BTCS";
  const std::string &approachString = this->approach_names[approach];
  std::string row = std::to_string(grid.getRow());
  std::string col = std::to_string(grid.getCol());
  std::string numIterations = std::to_string(iterations);

  std::string filename =
      approachString + "_" + row + "_" + col + "_" + numIterations + ".csv";

  while (std::filesystem::exists(filename)) {
    appendIdent += 1;
    appendIdentString = std::to_string(appendIdent);
    filename = approachString + "_" + row + "_" + col + "_" + numIterations +
               "-" + appendIdentString + ".csv";
  }

  file.open(filename);
  if (!file) {
    exit(1);
  }

  // adds lines at the beginning of verbose output csv that represent the
  // boundary conditions and their values -1 in case of closed boundary
  if (csv_output == CSV_OUTPUT_XTREME) {
    Eigen::IOFormat one_row(Eigen::StreamPrecision, Eigen::DontAlignCols, "",
                            " ");
    file << bc.getBoundarySideValues(BC_SIDE_LEFT).format(one_row)
         << std::endl; // boundary left
    file << bc.getBoundarySideValues(BC_SIDE_RIGHT).format(one_row)
         << std::endl; // boundary right
    file << bc.getBoundarySideValues(BC_SIDE_TOP).format(one_row)
         << std::endl; // boundary top
    file << bc.getBoundarySideValues(BC_SIDE_BOTTOM).format(one_row)
         << std::endl; // boundary bottom
    file << std::endl << std::endl;
  }

  file.close();

  return filename;
}

void Simulation::printConcentrationsCSV(const std::string &filename) {
  std::ofstream file;

  file.open(filename, std::ios_base::app);
  if (!file) {
    exit(1);
  }

  Eigen::IOFormat do_not_align(Eigen::StreamPrecision, Eigen::DontAlignCols);
  file << grid.getConcentrations().format(do_not_align) << std::endl;
  file << std::endl << std::endl;
  file.close();
}

void Simulation::run() {
  if (this->timestep == -1) {
    throw_invalid_argument("Timestep is not set!");
  }
  if (this->iterations == -1) {
    throw_invalid_argument("Number of iterations are not set!");
  }

  std::string filename;
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
      //     double percentage = (double)i / ((double)iterations *
      //     (double)innerIterations) * 100; if ((int)percentage % 10 == 0) {
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

    double beta = 0.5;

    // TODO this implementation is very inefficient!
    // a separate implementation that sets up a specific tridiagonal matrix for
    // Crank-Nicolson would be better
    Eigen::MatrixXd concentrations;
    Eigen::MatrixXd concentrationsFTCS;
    Eigen::MatrixXd concentrationsResult;
    for (int i = 0; i < iterations * innerIterations; i++) {
      if (console_output == CONSOLE_OUTPUT_VERBOSE && i > 0) {
        printConcentrationsConsole();
      }
      if (csv_output >= CSV_OUTPUT_VERBOSE) {
        printConcentrationsCSV(filename);
      }

      concentrations = grid.getConcentrations();
      FTCS(this->grid, this->bc, this->timestep, this->numThreads);
      concentrationsFTCS = grid.getConcentrations();
      grid.setConcentrations(concentrations);
      BTCS_Thomas(this->grid, this->bc, this->timestep, this->numThreads);
      concentrationsResult =
          beta * concentrationsFTCS + (1 - beta) * grid.getConcentrations();
      grid.setConcentrations(concentrationsResult);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  auto milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

  if (this->console_output > CONSOLE_OUTPUT_OFF) {
    printConcentrationsConsole();
  }
  if (this->csv_output > CSV_OUTPUT_OFF) {
    printConcentrationsCSV(filename);
  }
  if (this->time_measure > TIME_MEASURE_OFF) {
    const std::string &approachString = this->approach_names[approach];
    const std::string dimString = std::to_string(grid.getDim()) + "D";
    std::cout << approachString << dimString << ":: run() finished in "
              << milliseconds.count() << "ms" << std::endl;
  }
}
