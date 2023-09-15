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

template <class T> void Simulation<T>::setTimestep(T timestep) {
  if (timestep <= 0) {
    throw_invalid_argument("Timestep has to be greater than zero.");
  }

  if (approach == FTCS_APPROACH || approach == CRANK_NICOLSON_APPROACH) {

    const T deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();
    // will be 0 if 1D, else ...
    const T deltaRowSquare = grid.getDim() != 1
                                 ? grid.getDeltaRow() * grid.getDeltaRow()
                                 : deltaColSquare;
    const T minDeltaSquare =
        (deltaRowSquare < deltaColSquare) ? deltaRowSquare : deltaColSquare;

    T maxAlpha = std::numeric_limits<T>::quiet_NaN();

    // determine maximum alpha
    if (grid.getDim() == 2) {

      const T maxAlphaX = grid.getAlphaX().maxCoeff();
      const T maxAlphaY = grid.getAlphaY().maxCoeff();
      maxAlpha = (maxAlphaX > maxAlphaY) ? maxAlphaX : maxAlphaY;

    } else if (grid.getDim() == 1) {
      maxAlpha = grid.getAlpha().maxCoeff();

    } else {
      throw_invalid_argument("Critical error: Undefined number of dimensions!");
    }

    const std::string dim = std::to_string(grid.getDim()) + "D";

    // Courant-Friedrichs-Lewy condition
    T cfl = minDeltaSquare / (4 * maxAlpha);

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

template <class T> void Simulation<T>::run() {
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

    constexpr T beta = 0.5;

    // TODO this implementation is very inefficient!
    // a separate implementation that sets up a specific tridiagonal matrix for
    // Crank-Nicolson would be better
    Eigen::MatrixX<T> concentrations;
    Eigen::MatrixX<T> concentrationsFTCS;
    Eigen::MatrixX<T> concentrationsResult;
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

template void Simulation<double>::setTimestep(double timestep);
template void Simulation<float>::setTimestep(float timestep);

template void Simulation<double>::run();
template void Simulation<float>::run();
