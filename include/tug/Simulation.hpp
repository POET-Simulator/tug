/**
 * @file Simulation.hpp
 * @brief API of Simulation class, that holds all information regarding a
 * specific simulation run like its timestep, number of iterations and output
 * options. Simulation object also holds a predefined Grid and Boundary object.
 *
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "Boundary.hpp"
#include "Grid.hpp"
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs() 1
#endif

namespace tug {

/**
 * @brief Enum defining the two implemented solution approaches.
 *
 */
enum APPROACH {
  FTCS_APPROACH, // Forward Time-Centered Space
  BTCS_APPROACH, // Backward Time-Centered Space solved with EigenLU solver
  CRANK_NICOLSON_APPROACH
};

/**
 * @brief Enum defining the Linear Equation solvers
 *
 */
enum SOLVER {
  EIGEN_LU_SOLVER,        // EigenLU solver
  THOMAS_ALGORITHM_SOLVER // Thomas Algorithm solver; more efficient for
                          // tridiagonal matrices
};

/**
 * @brief Enum holding different options for .csv output.
 *
 */
enum CSV_OUTPUT {
  CSV_OUTPUT_OFF,     // do not produce csv output
  CSV_OUTPUT_ON,      // produce csv output with last concentration matrix
  CSV_OUTPUT_VERBOSE, // produce csv output with all concentration matrices
  CSV_OUTPUT_XTREME   // csv output like VERBOSE but additional boundary
                      // conditions at beginning
};

/**
 * @brief Enum holding different options for console output.
 *
 */
enum CONSOLE_OUTPUT {
  CONSOLE_OUTPUT_OFF,    // do not print any output to console
  CONSOLE_OUTPUT_ON,     // print before and after concentrations to console
  CONSOLE_OUTPUT_VERBOSE // print all concentration matrices to console
};

/**
 * @brief Enum holding different options for time measurement.
 *
 */
enum TIME_MEASURE {
  TIME_MEASURE_OFF, // do not print any time measures
  TIME_MEASURE_ON   // print time measure after last iteration
};

/**
 * @brief The class forms the interface for performing the diffusion simulations
 * and contains all the methods for controlling the desired parameters, such as
 * time step, number of simulations, etc.
 *
 */

template <class T> class Simulation {
public:
  /**
   * @brief Set up a simulation environment. The timestep and number of
   * iterations must be set. For the BTCS approach, the Thomas algorithm is used
   * as the default linear equation solver as this is faster for tridiagonal
   *        matrices. CSV output, console output and time measure are off by
   * default. Also, the number of cores is set to the maximum number of cores -1
   * by default.
   *
   * @param grid Valid grid object
   * @param bc Valid boundary condition object
   * @param approach Approach to solving the problem. Either FTCS or BTCS.
   */
  Simulation(Grid<T> &_grid, Boundary<T> &_bc, APPROACH _approach)
      : grid(_grid), bc(_bc), approach(_approach){};

  /**
   * @brief Set the option to output the results to a CSV file. Off by default.
   *
   *
   * @param csv_output Valid output option. The following options can be set
   *                   here:
   *                     - CSV_OUTPUT_OFF: do not produce csv output
   *                     - CSV_OUTPUT_ON: produce csv output with last
   *                       concentration matrix
   *                     - CSV_OUTPUT_VERBOSE: produce csv output with all
   *                       concentration matrices
   *                     - CSV_OUTPUT_XTREME: produce csv output with all
   *                       concentration matrices and simulation environment
   */
  void setOutputCSV(CSV_OUTPUT csv_output) {
    if (csv_output < CSV_OUTPUT_OFF && csv_output > CSV_OUTPUT_VERBOSE) {
      throw std::invalid_argument("Invalid CSV output option given!");
    }

    this->csv_output = csv_output;
  }

  /**
   * @brief Set the options for outputting information to the console. Off by
   * default.
   *
   * @param console_output Valid output option. The following options can be set
   *                       here:
   *                        - CONSOLE_OUTPUT_OFF: do not print any output to
   * console
   *                        - CONSOLE_OUTPUT_ON: print before and after
   * concentrations to console
   *                        - CONSOLE_OUTPUT_VERBOSE: print all concentration
   * matrices to console
   */
  void setOutputConsole(CONSOLE_OUTPUT console_output) {
    if (console_output < CONSOLE_OUTPUT_OFF &&
        console_output > CONSOLE_OUTPUT_VERBOSE) {
      throw std::invalid_argument("Invalid console output option given!");
    }

    this->console_output = console_output;
  }

  /**
   * @brief Set the Time Measure option. Off by default.
   *
   * @param time_measure The following options are allowed:
   *                     - TIME_MEASURE_OFF: Time of simulation is not printed
   * to console
   *                     - TIME_MEASURE_ON: Time of simulation run is printed to
   * console
   */
  void setTimeMeasure(TIME_MEASURE time_measure) {
    if (time_measure < TIME_MEASURE_OFF && time_measure > TIME_MEASURE_ON) {
      throw std::invalid_argument("Invalid time measure option given!");
    }

    this->time_measure = time_measure;
  }

  /**
   * @brief Setting the time step for each iteration step. Time step must be
   *        greater than zero. Setting the timestep is required.
   *
   * @param timestep Valid timestep greater than zero.
   */
  void setTimestep(T timestep);

  /**
   * @brief Currently set time step is returned.
   *
   * @return double timestep
   */
  T getTimestep() const { return this->timestep; }

  /**
   * @brief Set the desired iterations to be calculated. A value greater
   *        than zero must be specified here. Setting iterations is required.
   *
   * @param iterations Number of iterations to be simulated.
   */
  void setIterations(int iterations) {
    if (iterations <= 0) {
      throw std::invalid_argument(
          "Number of iterations must be greater than zero.");
    }
    this->iterations = iterations;
  }

  /**
   * @brief Set the desired linear equation solver to be used for BTCS approach.
   * Without effect in case of FTCS approach.
   *
   * @param solver Solver to be used. Default is Thomas Algorithm as it is more
   * efficient for tridiagonal Matrices.
   */
  void setSolver(SOLVER solver) {
    if (this->approach == FTCS_APPROACH) {
      std::cerr
          << "Warning: Solver was set, but FTCS approach initialized. Setting "
             "the solver "
             "is thus without effect."
          << std::endl;
    }

    this->solver = solver;
  }

  /**
   * @brief Set the number of desired openMP Threads.
   *
   * @param num_threads Number of desired threads. Must have a value between
   *                    1 and the maximum available number of processors. The
   * maximum number of processors is set as the default case during Simulation
   * construction.
   */
  void setNumberThreads(int num_threads) {
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
   * @brief Return the currently set iterations to be calculated.
   *
   * @return int Number of iterations.
   */
  int getIterations() const { return this->iterations; }

  /**
   * @brief Outputs the current concentrations of the grid on the console.
   *
   */
  inline void printConcentrationsConsole() const {
    std::cout << grid.getConcentrations() << std::endl;
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
    file << grid.getConcentrations().format(do_not_align) << std::endl;
    file << std::endl << std::endl;
    file.close();
  }

  /**
   * @brief Method starts the simulation process with the previously set
   *        parameters.
   */
  void run();

private:
  T timestep{-1};
  int iterations{-1};
  int innerIterations{1};
  int numThreads{omp_get_num_procs()};
  CSV_OUTPUT csv_output{CSV_OUTPUT_OFF};
  CONSOLE_OUTPUT console_output{CONSOLE_OUTPUT_OFF};
  TIME_MEASURE time_measure{TIME_MEASURE_OFF};

  Grid<T> &grid;
  Boundary<T> &bc;
  APPROACH approach;
  SOLVER solver{THOMAS_ALGORITHM_SOLVER};

  const std::vector<std::string> approach_names = {"FTCS", "BTCS", "CRNI"};
};
} // namespace tug
#endif // SIMULATION_H_
