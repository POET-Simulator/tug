#pragma once

#include "tug/Boundary.hpp"
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

/**
 * @brief A base class for simulation grids.
 *
 * This class provides a base implementation for simulation grids, including
 * methods for setting and getting grid dimensions, domain sizes, and output
 * options. It also includes methods for running simulations, which must be
 * implemented by derived classes.
 *
 * @tparam T The type of the elements in the grid.
 */
template <typename T> class BaseSimulationGrid {
private:
  CSV_OUTPUT csv_output{CSV_OUTPUT::OFF};
  CONSOLE_OUTPUT console_output{CONSOLE_OUTPUT::OFF};
  TIME_MEASURE time_measure{TIME_MEASURE::OFF};

  int iterations{1};
  RowMajMatMap<T> concentrationMatrix;
  Boundary<T> boundaryConditions;

  const std::uint8_t dim;

  T delta_col;
  T delta_row;

protected:
  void applyInnerBoundaries() {
    const auto &inner_bc = boundaryConditions.getInnerBoundaries();
    if (inner_bc.empty()) {
      return;
    }

    for (const auto &[rowcol, value] : inner_bc) {
      concentrationMatrix(rowcol.first, rowcol.second) = value;
    }
  }

public:
  /**
   * @brief Constructs a BaseSimulationGrid from a given RowMajMat object.
   *
   * This constructor initializes a BaseSimulationGrid using the data, number of
   * rows, and number of columns from the provided RowMajMat object.
   *
   * @tparam T The type of the elements in the RowMajMat.
   * @param origin The RowMajMat object from which to initialize the
   * BaseSimulationGrid.
   */
  BaseSimulationGrid(RowMajMat<T> &origin)
      : BaseSimulationGrid(origin.data(), origin.rows(), origin.cols()) {}

  /**
   * @brief Constructs a BaseSimulationGrid object.
   *
   * @tparam T The type of the data elements.
   * @param data Pointer to the data array.
   * @param rows Number of rows in the grid.
   * @param cols Number of columns in the grid.
   *
   * Initializes the concentration_matrix with the provided data, rows, and
   * columns. Sets delta_col and delta_row to 1. Determines the dimension (dim)
   * based on the number of rows: if rows == 1, dim is set to 1; otherwise, dim
   * is set to 2.
   */
  BaseSimulationGrid(T *data, std::size_t rows, std::size_t cols)
      : concentrationMatrix(data, rows, cols), boundaryConditions(rows, cols),
        delta_col(1), delta_row(1), dim(rows == 1 ? 1 : 2) {}

  /**
   * @brief Constructs a BaseSimulationGrid with a single dimension.
   *
   * This constructor initializes a BaseSimulationGrid object with the provided
   * data and length. It assumes the grid has only one dimension.
   *
   * @param data Pointer to the data array.
   * @param length The length of the data array.
   */
  BaseSimulationGrid(T *data, std::size_t length)
      : BaseSimulationGrid(data, 1, length) {}

  /**
   * @brief Overloaded function call operator to access elements in a
   * one-dimensional grid.
   *
   * This operator provides access to elements in the concentration matrix using
   * a single index. It asserts that the grid is one-dimensional before
   * accessing the element.
   *
   * @tparam T The type of elements in the concentration matrix.
   * @param index The index of the element to access.
   * @return A reference to the element at the specified index in the
   * concentration matrix.
   */
  constexpr T &operator()(std::size_t index) {
    tug_assert(dim == 1, "Grid is not one dimensional, use 2D index operator!");

    return concentrationMatrix(index);
  }

  /**
   * @brief Overloaded function call operator to access elements in a 2D
   * concentration matrix.
   *
   * This operator allows accessing elements in the concentration matrix using
   * row and column indices. It asserts that the grid is two-dimensional before
   * accessing the element.
   *
   * @param row The row index of the element to access.
   * @param col The column index of the element to access.
   * @return A reference to the element at the specified row and column in the
   * concentration matrix.
   */
  constexpr T &operator()(std::size_t row, std::size_t col) {
    tug_assert(dim == 2, "Grid is not two dimensional, use 1D index operator!");

    return concentrationMatrix(row, col);
  }

  /**
   * @brief Retrieves the concentration matrix.
   *
   * @tparam T The data type of the elements in the concentration matrix.
   * @return RowMajMat<T>& Reference to the concentration matrix.
   */
  RowMajMatMap<T> &getConcentrationMatrix() { return concentrationMatrix; }

  const RowMajMatMap<T> &getConcentrationMatrix() const {
    return concentrationMatrix;
  }

  /**
   * @brief Retrieves the boundary conditions for the simulation.
   *
   * @tparam T The type parameter for the Boundary class.
   * @return Boundary<T>& A reference to the boundary conditions.
   */
  Boundary<T> &getBoundaryConditions() { return boundaryConditions; }

  const Boundary<T> &getBoundaryConditions() const {
    return boundaryConditions;
  }

  /**
   * @brief Retrieves the dimension value.
   *
   * @return The dimension value as an 8-bit unsigned integer.
   */
  std::uint8_t getDim() const { return dim; }

  /**
   * @brief Returns the number of rows in the concentration matrix.
   *
   * @return std::size_t The number of rows in the concentration matrix.
   */
  std::size_t rows() const { return concentrationMatrix.rows(); }

  /**
   * @brief Get the number of columns in the concentration matrix.
   *
   * @return std::size_t The number of columns in the concentration matrix.
   */
  std::size_t cols() const { return concentrationMatrix.cols(); }

  /**
   * @brief Returns the cell size in meter of the x-direction.
   *
   * This function returns the value of the delta column, which is used
   * to represent the difference or change in the column value.
   *
   * @return T The cell size in meter of the x-direction.
   */
  T deltaCol() const { return delta_col; }

  /**
   * @brief Returns the cell size in meter of the y-direction.
   *
   * This function asserts that the grid is two-dimensional. If the grid is not
   * two-dimensional, an assertion error is raised with the message "Grid is not
   * two dimensional, there is no delta in y-direction!".
   *
   * @return The cell size in meter of the y-direction.
   */
  T deltaRow() const {
    tug_assert(
        dim == 2,
        "Grid is not two dimensional, there is no delta in y-direction!");

    return delta_row;
  }

  /**
   * @brief Computes the domain size in the X direction.
   *
   * This function calculates the size of the domain in the X direction by
   * multiplying the column spacing (delta_col) by the number of columns (cols).
   *
   * @return The size of the domain in the X direction.
   */
  T domainX() const { return delta_col * cols(); }

  /**
   * @brief Returns the size of the domain in the y-direction.
   *
   * This function calculates the size of the domain in the y-direction
   * by multiplying the row spacing (delta_row) by the number of rows.
   * It asserts that the grid is two-dimensional before performing the
   * calculation.
   *
   * @return The size of the domain in the y-direction.
   */
  T domainY() const {
    tug_assert(
        dim == 2,
        "Grid is not two dimensional, there is no domain in y-direction!");

    return delta_row * rows();
  }

  /**
   * @brief Sets the domain length for a one-dimensional grid.
   *
   * This function sets the domain length for a one-dimensional grid and
   * calculates the column width (delta_col) based on the given domain length
   * and the number of columns. It asserts that the grid is one-dimensional and
   * that the given domain length is positive.
   *
   * @param domain_length The length of the domain. Must be positive.
   */
  void setDomain(T domain_length) {
    tug_assert(dim == 1, "Grid is not one dimensional, use 2D domain setter!");
    tug_assert(domain_length > 0, "Given domain length is not positive!");

    delta_col = domain_length / cols();
  }

  /**
   * @brief Sets the domain size for a 2D grid simulation.
   *
   * This function sets the domain size in the x and y directions for a
   * two-dimensional grid simulation. It asserts that the grid is indeed
   * two-dimensional and that the provided domain sizes are positive.
   *
   * @tparam T The type of the domain size parameters.
   * @param domain_row The size of the domain in the y-direction.
   * @param domain_col The size of the domain in the x-direction.
   */
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
   * @brief Retrieves the CSV output.
   *
   * This function returns the CSV output associated with the simulation.
   *
   * @return CSV_OUTPUT The CSV output of the simulation.
   */
  constexpr CSV_OUTPUT getOutputCSV() const { return this->csv_output; }

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
   * @brief Retrieves the console output.
   *
   * This function returns the current state of the console output.
   *
   * @return CONSOLE_OUTPUT The current console output.
   */
  constexpr CONSOLE_OUTPUT getOutputConsole() const {
    return this->console_output;
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
   * @brief Retrieves the current time measurement.
   *
   * @return TIME_MEASURE The current time measurement.
   */
  constexpr TIME_MEASURE getTimeMeasure() const { return this->time_measure; }

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

  virtual void setTimestep(T timestep) = 0;
};
} // namespace tug
