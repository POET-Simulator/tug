/**
 * @file Diffusion.hpp
 * @brief API of Diffusion class, that holds all information regarding a
 * specific simulation run like its timestep, number of iterations and output
 * options. Diffusion object also holds a predefined Grid and Boundary object.
 *
 */

#pragma once

#include "tug/Core/Matrix.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <tug/Core/BaseSimulation.hpp>
#include <tug/Core/Numeric/BTCS.hpp>
#include <tug/Core/Numeric/FTCS.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs() 1
#endif

namespace tug {

/**
 * @brief Enum defining the implemented solution approaches.
 *
 */
enum APPROACH {
  FTCS_APPROACH,          /*!< Forward Time-Centered Space */
  BTCS_APPROACH,          /*!< Backward Time-Centered Space */
  CRANK_NICOLSON_APPROACH /*!< Crank-Nicolson method */
};

/**
 * @brief Enum defining the Linear Equation solvers
 *
 */
enum SOLVER {
  EIGEN_LU_SOLVER,        /*!<  EigenLU solver */
  THOMAS_ALGORITHM_SOLVER /*!< Thomas Algorithm solver; more efficient for
                             tridiagonal matrices */
};

/**
 * @brief The class forms the interface for performing the diffusion simulations
 * and contains all the methods for controlling the desired parameters, such as
 * time step, number of simulations, etc.
 *
 * @tparam T the type of the internal data structures for grid, boundary
 * condition and timestep
 * @tparam approach Set the SLE scheme to be used
 * @tparam solver Set the solver to be used
 */
template <class T, APPROACH approach = BTCS_APPROACH,
          SOLVER solver = THOMAS_ALGORITHM_SOLVER>
class Diffusion : public BaseSimulationGrid<T> {
private:
  T timestep{-1};
  int innerIterations{1};
  int numThreads{omp_get_num_procs()};

  RowMajMat<T> alphaX;
  RowMajMat<T> alphaY;

  const std::vector<std::string> approach_names = {"FTCS", "BTCS", "CRNI"};

  static constexpr T DEFAULT_ALPHA = 1E-8;

  void init_alpha() {
    this->alphaX =
        RowMajMat<T>::Constant(this->rows(), this->cols(), DEFAULT_ALPHA);
    if (this->getDim() == 2) {
      this->alphaY =
          RowMajMat<T>::Constant(this->rows(), this->cols(), DEFAULT_ALPHA);
    }
  }

public:
  Diffusion(RowMajMat<T> &origin) : BaseSimulationGrid<T>(origin) {
    init_alpha();
  }

  Diffusion(T *data, int rows, int cols)
      : BaseSimulationGrid<T>(data, rows, cols) {
    init_alpha();
  }

  Diffusion(T *data, std::size_t length) : BaseSimulationGrid<T>(data, length) {
    init_alpha();
  }

  RowMajMat<T> &getAlphaX() { return alphaX; }

  RowMajMat<T> &getAlphaY() {
    tug_assert(
        this->getDim(),
        "Grid is not two dimensional, there is no domain in y-direction!");

    return alphaY;
  }

  void setAlphaX(const RowMajMat<T> &alphaX) { this->alphaX = alphaX; }

  void setAlphaY(const RowMajMat<T> &alphaY) {
    tug_assert(
        this->getDim(),
        "Grid is not two dimensional, there is no domain in y-direction!");

    this->alphaY = alphaY;
  }

  // /**
  //  * @brief Set up a simulation environment. The timestep and number of
  //  * iterations must be set. For the BTCS approach, the Thomas algorithm is
  //  used
  //  * as the default linear equation solver as this is faster for tridiagonal
  //  *        matrices. CSV output, console output and time measure are off by
  //  * default. Also, the number of cores is set to the maximum number of cores
  //  -1
  //  * by default.
  //  *
  //  * @param grid Valid grid object
  //  * @param bc Valid boundary condition object
  //  * @param approach Approach to solving the problem. Either FTCS or BTCS.
  //  */
  // Diffusion(Grid<T> &_grid, Boundary<T> &_bc) : grid(_grid), bc(_bc) {};

  /**
   * @brief Setting the time step for each iteration step. Time step must be
   *        greater than zero. Setting the timestep is required.
   *
   * @param timestep Valid timestep greater than zero.
   */
  void setTimestep(T timestep) override {
    tug_assert(timestep > 0, "Timestep has to be greater than zero.");

    if constexpr (approach == FTCS_APPROACH ||
                  approach == CRANK_NICOLSON_APPROACH) {
      T cfl;
      if (this->getDim() == 1) {

        const T deltaSquare = this->deltaCol();
        const T maxAlpha = this->alphaX.maxCoeff();

        // Courant-Friedrichs-Lewy condition
        cfl = deltaSquare / (4 * maxAlpha);
      } else if (this->getDim() == 2) {
        const T deltaColSquare = this->deltaCol() * this->deltaCol();
        // will be 0 if 1D, else ...
        const T deltaRowSquare = this->deltaRow() * this->deltaRow();
        const T minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);

        const T maxAlpha =
            std::max(this->alphaX.maxCoeff(), this->alphaY.maxCoeff());

        cfl = minDeltaSquare / (4 * maxAlpha);
      }
      const std::string dim = std::to_string(this->getDim()) + "D";

      const std::string &approachPrefix = this->approach_names[approach];
      std::cout << approachPrefix << "_" << dim << " :: CFL condition: " << cfl
                << std::endl;
      std::cout << approachPrefix << "_" << dim
                << " :: required dt=" << timestep << std::endl;

      if (timestep > cfl) {

        this->innerIterations = (int)ceil(timestep / cfl);
        this->timestep = timestep / (double)innerIterations;

        std::cerr << "Warning :: Timestep was adjusted, because of stability "
                     "conditions. Time duration was approximately preserved by "
                     "adjusting internal number of iterations."
                  << std::endl;
        std::cout << approachPrefix << "_" << dim << " :: Required "
                  << this->innerIterations
                  << " inner iterations with dt=" << this->timestep
                  << std::endl;

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

  /**
   * @brief Currently set time step is returned.
   *
   * @return double timestep
   */
  T getTimestep() const { return this->timestep; }

  /**
   * @brief Set the number of desired openMP Threads.
   *
   * @param num_threads Number of desired threads. Must have a value between
   *                    1 and the maximum available number of processors. The
   * maximum number of processors is set as the default case during Simulation
   * construction.
   */
  void setNumberThreads(int numThreads) {
    if (numThreads > 0 && numThreads <= omp_get_num_procs()) {
      this->numThreads = numThreads;
    } else {
      int maxThreadNumber = omp_get_num_procs();
      throw std::invalid_argument(
          "Number of threads exceeds the number of processor cores (" +
          std::to_string(maxThreadNumber) + ") or is less than 1.");
    }
  }

  /**
   * @brief Outputs the current concentrations of the grid on the console.
   *
   */
  void printConcentrationsConsole() const {
    std::cout << this->getConcentrationMatrix() << std::endl;
    std::cout << std::endl;
  }

  /**
   * @brief Creates a CSV file with a name containing the current simulation
   *        parameters. If the data name already exists, an additional counter
   * is appended to the name. The name of the file is built up as follows:
   *        <approach> + <number rows> + <number columns> + <number of
   * iterations>+<counter>.csv
   *
   * @return string Filename with configured simulation parameters.
   */
  std::string createCSVfile() const {
    std::ofstream file;
    int appendIdent = 0;
    std::string appendIdentString;

    // string approachString = (approach == 0) ? "FTCS" : "BTCS";
    const std::string &approachString = this->approach_names[approach];
    std::string row = std::to_string(this->rows());
    std::string col = std::to_string(this->cols());
    std::string numIterations = std::to_string(this->getIterations());

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
    if (this->getOutputCSV() == CSV_OUTPUT::XTREME) {
      const auto &bc = this->getBoundaryConditions();

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

  /**
   * @brief Writes the currently calculated concentration values of the grid
   *        into the CSV file with the passed filename.
   *
   * @param filename Name of the file to which the concentration values are
   *                 to be written.
   */
  void printConcentrationsCSV(const std::string &filename) const {
    std::ofstream file;

    file.open(filename, std::ios_base::app);
    if (!file) {
      exit(1);
    }

    Eigen::IOFormat do_not_align(Eigen::StreamPrecision, Eigen::DontAlignCols);
    file << this->getConcentrationMatrix().format(do_not_align) << std::endl;
    file << std::endl << std::endl;
    file.close();
  }

  /**
   * @brief Method starts the simulation process with the previously set
   *        parameters.
   */
  void run() override {
    tug_assert(this->getTimestep() > 0, "Timestep is not set!");
    tug_assert(this->getIterations() > 0, "Number of iterations are not set!");

    std::string filename;
    if (this->getOutputConsole() > CONSOLE_OUTPUT::OFF) {
      printConcentrationsConsole();
    }
    if (this->getOutputCSV() > CSV_OUTPUT::OFF) {
      filename = createCSVfile();
    }

    auto begin = std::chrono::high_resolution_clock::now();

    SimulationInput<T> sim_input = {.concentrations =
                                        this->getConcentrationMatrix(),
                                    .alphaX = this->getAlphaX(),
                                    .alphaY = this->getAlphaY(),
                                    .boundaries = this->getBoundaryConditions(),
                                    .dim = this->getDim(),
                                    .timestep = this->getTimestep(),
                                    .rowMax = this->rows(),
                                    .colMax = this->cols(),
                                    .deltaRow = this->deltaRow(),
                                    .deltaCol = this->deltaCol()};

    if constexpr (approach == FTCS_APPROACH) { // FTCS case
      for (int i = 0; i < this->getIterations() * innerIterations; i++) {
        if (this->getOutputConsole() == CONSOLE_OUTPUT::VERBOSE && i > 0) {
          printConcentrationsConsole();
        }
        if (this->getOutputCSV() >= CSV_OUTPUT::VERBOSE) {
          printConcentrationsCSV(filename);
        }

        FTCS(sim_input, this->numThreads);

        // if (i % (iterations * innerIterations / 100) == 0) {
        //     double percentage = (double)i / ((double)iterations *
        //     (double)innerIterations) * 100; if ((int)percentage % 10 == 0) {
        //         cout << "Progress: " << percentage << "%" << endl;
        //     }
        // }
      }

    } else if constexpr (approach == BTCS_APPROACH) { // BTCS case

      if constexpr (solver == EIGEN_LU_SOLVER) {
        for (int i = 0; i < this->getIterations(); i++) {
          if (this->getOutputConsole() == CONSOLE_OUTPUT::VERBOSE && i > 0) {
            printConcentrationsConsole();
          }
          if (this->getOutputCSV() >= CSV_OUTPUT::VERBOSE) {
            printConcentrationsCSV(filename);
          }

          BTCS_LU(sim_input, this->numThreads);
        }
      } else if constexpr (solver == THOMAS_ALGORITHM_SOLVER) {
        for (int i = 0; i < this->getIterations(); i++) {
          if (this->getOutputConsole() == CONSOLE_OUTPUT::VERBOSE && i > 0) {
            printConcentrationsConsole();
          }
          if (this->getOutputCSV() >= CSV_OUTPUT::VERBOSE) {
            printConcentrationsCSV(filename);
          }

          BTCS_Thomas(sim_input, this->numThreads);
        }
      }
    } else if constexpr (approach ==
                         CRANK_NICOLSON_APPROACH) { // Crank-Nicolson case

      constexpr T beta = 0.5;

      // TODO this implementation is very inefficient!
      // a separate implementation that sets up a specific tridiagonal matrix
      // for Crank-Nicolson would be better
      RowMajMat<T> concentrations;
      RowMajMat<T> concentrationsFTCS;
      RowMajMat<T> concentrationsResult;
      for (int i = 0; i < this->getIterations() * innerIterations; i++) {
        if (this->getOutputConsole() == CONSOLE_OUTPUT::VERBOSE && i > 0) {
          printConcentrationsConsole();
        }
        if (this->getOutputCSV() >= CSV_OUTPUT::VERBOSE) {
          printConcentrationsCSV(filename);
        }

        concentrations = this->getConcentrationMatrix();
        FTCS(this->grid, this->bc, this->timestep, this->numThreads);
        concentrationsFTCS = this->getConcentrationMatrix();
        this->getConcentrationMatrix() = concentrations;
        BTCS_Thomas(sim_input, this->numThreads);
        concentrationsResult = beta * concentrationsFTCS +
                               (1 - beta) * this->getConcentrationMatrix();
        this->getConcentrationMatrix() = concentrationsResult;
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    if (this->getOutputConsole() > CONSOLE_OUTPUT::OFF) {
      printConcentrationsConsole();
    }
    if (this->getOutputCSV() > CSV_OUTPUT::OFF) {
      printConcentrationsCSV(filename);
    }
    if (this->getTimeMeasure() > TIME_MEASURE::OFF) {
      const std::string &approachString = this->approach_names[approach];
      const std::string dimString = std::to_string(this->getDim()) + "D";
      std::cout << approachString << dimString << ":: run() finished in "
                << milliseconds.count() << "ms" << std::endl;
    }
  }
};
} // namespace tug
