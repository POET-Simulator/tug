#pragma once

#include <stdexcept>

namespace tug {

/**
 * @brief Enum holding different options for .csv output.
 *
 */
enum CSV_OUTPUT {
  CSV_OUTPUT_OFF,     /*!< do not produce csv output */
  CSV_OUTPUT_ON,      /*!< produce csv output with last concentration matrix */
  CSV_OUTPUT_VERBOSE, /*!< produce csv output with all concentration matrices */
  CSV_OUTPUT_XTREME   /*!< csv output like VERBOSE but additional boundary
                         conditions at beginning */
};

/**
 * @brief Enum holding different options for console output.
 *
 */
enum CONSOLE_OUTPUT {
  CONSOLE_OUTPUT_OFF, /*!< do not print any output to console */
  CONSOLE_OUTPUT_ON,  /*!< print before and after concentrations to console */
  CONSOLE_OUTPUT_VERBOSE /*!< print all concentration matrices to console */
};

/**
 * @brief Enum holding different options for time measurement.
 *
 */
enum TIME_MEASURE {
  TIME_MEASURE_OFF, /*!< do not print any time measures */
  TIME_MEASURE_ON   /*!< print time measure after last iteration */
};

class BaseSimulation {
protected:
  CSV_OUTPUT csv_output{CSV_OUTPUT_OFF};
  CONSOLE_OUTPUT console_output{CONSOLE_OUTPUT_OFF};
  TIME_MEASURE time_measure{TIME_MEASURE_OFF};

  int iterations{1};

public:
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