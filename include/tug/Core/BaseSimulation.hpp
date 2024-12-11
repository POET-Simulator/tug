#pragma once

#include <cstddef>
#include <cstdint>
#include <tug/Core/Matrix.hpp>
#include <tug/Core/TugUtils.hpp>

namespace tug {

/**
 * @brief Enum holding different options for .csv output.
 *
 */
enum class CSV_OUTPUT {
  OFF,     /*!< do not produce csv output */
  ON,      /*!< produce csv output with last concentration matrix */
  VERBOSE, /*!< produce csv output with all concentration matrices */
  XTREME   /*!< csv output like VERBOSE but additional boundary
                         conditions at beginning */
};

/**
 * @brief Enum holding different options for console output.
 *
 */
enum class CONSOLE_OUTPUT {
  OFF,    /*!< do not print any output to console */
  ON,     /*!< print before and after concentrations to console */
  VERBOSE /*!< print all concentration matrices to console */
};

/**
 * @brief Enum holding different options for time measurement.
 *
 */
enum class TIME_MEASURE {
  OFF, /*!< do not print any time measures */
  ON   /*!< print time measure after last iteration */
};

template <typename T> class BaseSimulationGrid {
protected:
  CSV_OUTPUT csv_output{CSV_OUTPUT::OFF};
  CONSOLE_OUTPUT console_output{CONSOLE_OUTPUT::OFF};
  TIME_MEASURE time_measure{TIME_MEASURE::OFF};

  int iterations{1};
  RowMajMatMap<T> concentration_matrix;

  const std::uint8_t dim;

  T delta_col;
  T delta_row;

public:
  BaseSimulationGrid(T *data, std::size_t length)
      : BaseSimulationGrid(data, 1, length) {}

  template <typename EigenType>
  BaseSimulationGrid(const EigenType &origin)
      : BaseSimulationGrid(origin.data(), origin.rows(), origin.cols()) {}

  BaseSimulationGrid(T *data, std::size_t rows, std::size_t cols)
      : concentration_matrix(data, rows, cols), delta_col(1), delta_row(1),
        dim(rows == 1 ? 1 : 2) {}

  std::size_t rows() const { return concentration_matrix.rows(); }
  std::size_t cols() const { return concentration_matrix.cols(); }

  T deltaCol() const { return delta_col; }
  T deltaRow() const {
    tug_assert(
        dim == 1,
        "Grid is not two dimensional, there is no delta in y-direction!");

    return delta_row;
  }

  void setDomain(T domain_length) {
    tug_assert(dim == 1, "Grid is not one dimensional, use 2D domain setter!");
    tug_assert(domain_length > 0, "Given domain length is not positive!");

    delta_col = domain_length / cols();
  }

  void setDomain(T domain_row, T domain_col) {
    tug_assert(dim == 2, "Grid is not two dimensional, use 1D domain setter!");
    tug_assert(domain_col > 0,
               "Given domain size in x-direction is not positive!");
    tug_assert(domain_row > 0,
               "Given domain size in y-direction is not positive!");

    delta_row = domain_row / rows();
    delta_col = domain_col / cols();
  }

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
  void setOutputCSV(CSV_OUTPUT csv_output) { this->csv_output = csv_output; }

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
    this->time_measure = time_measure;
  }

  /**
   * @brief Set the desired iterations to be calculated. A value greater
   *        than zero must be specified here. Setting iterations is required.
   *
   * @param iterations Number of iterations to be simulated.
   */
  void setIterations(int iterations) {
    tug_assert(iterations > 0,
               "Number of iterations must be greater than zero.");

    this->iterations = iterations;
  }

  /**
   * @brief Return the currently set iterations to be calculated.
   *
   * @return int Number of iterations.
   */
  int getIterations() const { return this->iterations; }

  /**
   * @brief Method starts the simulation process with the previously set
   *        parameters.
   */
  virtual void run() = 0;
};
} // namespace tug