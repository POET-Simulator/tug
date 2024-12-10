/**
 * @file Advection.hpp
 * @brief API of Advection class, holding information for a simulation of
 * advection. Holds a predifined Grid object, Boundary object and Velocities
 * object
 */

#pragma once

#include "tug/Core/Matrix.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <tug/Boundary.hpp>
#include <tug/Grid.hpp>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <tug/Core/Numeric/BTCS.hpp>
#include <tug/Core/Numeric/FTCS.hpp>
#include <tug/Core/TugUtils.hpp>
#include <tug/Diffusion/Diffusion.hpp>

#include <tug/Advection/Velocities.hpp>

using namespace Eigen;
namespace tug {
template <class T> class Advection {
public:
  /**
   * @brief Construct a new Advection object, used to calculate material
   * transport. A timestep and number of iterations must be set. A transient
   * case can be selected by initializing Steady=false. With each timestep the
   * Velocities object will also be updated.
   * A steady case can be selected by initializing Steady=true. The velocities
   * object will not be updated. Velocities can be simulated to convergence
   * beforehand. Porosity can be set, the default is 1. CSV Output is off by
   * default.
   *
   * @param grid Valid grid object
   * @param bc Valid Boundary object
   * @param Steady Used to choose between Steady and Transient case. Either true
   * or false
   */
  Advection(Velocities<T> &_velocities, Grid<T> &_grid, Boundary<T> &_bc,
            bool Steady)
      : velocities(_velocities), grid(_grid), bc(_bc),
        outx(_velocities.getOutx()), outy(_velocities.getOuty()),
        Steady(Steady) {};

  /**
   * @brief Sets the porosity of the medium
   *
   * @param porosity new porosity value
   */
  void setPorosity(T porosity) {
    if (porosity < 0 || porosity > 1) {
      throw std::invalid_argument(
          "Porosity must be a value between 0 and 1 (inclusive)");
    }
    this->porosity = porosity;
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
          "Number of iterations must be greater than zero. Provided value: " +
          std::to_string(iterations));
    }
    this->iterations = iterations;
  };

  /**
   * @brief Set the size of the timestep. Must be greater than zero
   *
   * @param timestep Size of the timestep
   */
  void setTimestep(T timestep) {
    if (timestep <= 0) {
      throw std::invalid_argument(
          "Timestep must be greater than zero. Provided value: " +
          std::to_string(timestep));
    } else {
      this->timestep = timestep;
    }
  }

  /**
   * @brief Set the number of desired openMP Threads.
   *
   * @param num_threads Number of desired threads. Must have a value between
   *                    1 and the maximum available number of processors. The
   * maximum number of processors is set as the default case during Advection
   * construction.
   */
  void setNumberThreads(int num_threads) {
    if (num_threads > 0 && num_threads <= omp_get_num_procs()) {
      this->numThreads = num_threads;
    } else {
      int maxThreadNumber = omp_get_num_procs();
      if (num_threads > maxThreadNumber) {
        throw std::invalid_argument(
            "Number of threads exceeds the number of processor cores (" +
            std::to_string(maxThreadNumber) + ").");
      } else {
        throw std::invalid_argument("Number of threads is less than 1.");
      }
    }
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
    if (csv_output < CSV_OUTPUT_OFF || csv_output > CSV_OUTPUT_VERBOSE) {
      throw std::invalid_argument("Invalid CSV output option given: " +
  std::to_string(csv_output));
    }
  void setOutputCSV(CSV_OUTPUT csv_output) {
    if (csv_output < CSV_OUTPUT_OFF && csv_output > CSV_OUTPUT_VERBOSE) {
      throw std::invalid_argument("Invalid CSV output option given!");
    }

    this->csv_output = csv_output;
  }

  /**
   * @brief Creates a CSV file with a name containing the current simulation
   *        parameters. If the data name already exists, an additional counter
   * is appended to the name. The name of the file is built up as follows:
   *        <Information contained in file> + <number rows> + <number columns> +
   * <number of iterations>+<counter>.csv
   *
   * @return string Filename with configured simulation parameters.
   */
  std::string createCSVfile(std::string Type) const {
    std::ofstream file;
    int appendIdent = 0;
    std::string appendIdentString;

    std::string row = std::to_string(grid.getRow());
    std::string col = std::to_string(grid.getCol());
    std::string numIterations = std::to_string(iterations);

    std::string filename =
        Type + "_" + row + "_" + col + "_" + numIterations + ".csv";

    while (std::filesystem::exists(filename)) {
      appendIdent += 1;
      appendIdentString = std::to_string(appendIdent);
      filename = Type + "_" + row + "_" + col + "_" + numIterations + "-" +
                 appendIdentString + ".csv";
    }

    file.open(filename);
    if (!file) {
      throw std::runtime_error("Failed to open file: " + filename);
    }

    file.close();

    return filename;
  }
  /**
   * @brief Writes the currently calculated Concentration values of the grid
   *        into the CSV file with the passed filename.
   *
   * @param filename Name of the file to which the Concentration values are
   *                 to be written.
   */
  void printConcentrationCSV(const std::string &filename) const {
    std::ofstream file;

    file.open(filename, std::ios_base::app);
    if (!file) {
      throw std::runtime_error("Failed to open file: " + filename);
    }

    Eigen::IOFormat do_not_align(Eigen::StreamPrecision, Eigen::DontAlignCols);
    file << grid.getConcentrations().format(do_not_align) << std::endl;
    file << std::endl << std::endl;
    file.close();
  }

  /**
   * @brief function calculating material transport for one timestep
   */
  void adv() {
    int rows = grid.getRow();
    int cols = grid.getCol();
    T volume = grid.getDeltaRow() * grid.getDeltaCol();

    RowMajMat<T> &newConcentrations = grid.getConcentrations();

    // Calculate Courant-Levy-Frederich condition
    T maxFx = std::max(abs(outx.maxCoeff()), abs(outx.minCoeff()));
    T maxFy = std::max(abs(outy.maxCoeff()), abs(outy.minCoeff()));
    T maxF = std::max(maxFx, maxFy);

    if (maxF == 0) {
      throw std::runtime_error("Division by zero: maxF is zero.");
    }

    T cvf = abs((volume * porosity) / maxF);
    int innerSteps = (int)ceil(timestep / cvf);
    T innerTimestep = timestep / innerSteps;

    for (int k = 0; k < innerSteps; k++) {
      const Eigen::MatrixX<T> oldConcentrations = newConcentrations;
// Calculate sum of incoming/outgoing Flow*Concentration in x-direction in each
// cell
#pragma omp parallel for num_threads(numThreads) schedule(static)
      for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols + 1; j++) {
          if (j == 0) {
            if (bc.getBoundaryElementType(BC_SIDE_LEFT, i) != BC_TYPE_CLOSED) {
              if (outx(i, j) > 0) {
                // outx positive -> flow from border to cell i,j
                newConcentrations(i, j) +=
                    outx(i, j) * bc.getBoundaryElementValue(BC_SIDE_LEFT, i);
              } else if (outx(i, j) < 0) {
                // outx negative -> flow from i,j towards border
                newConcentrations(i, j) += outx(i, j) * oldConcentrations(i, j);
              }
            }
          } else if (j == cols) {
            if (bc.getBoundaryElementType(BC_SIDE_RIGHT, i) != BC_TYPE_CLOSED) {
              if (outx(i, j) > 0) {
                // outx positive-> flow from i,j-1 towards border
                newConcentrations(i, j - 1) -=
                    outx(i, j) * oldConcentrations(i, j - 1);
              } else if (outx(i, j) < 0) {
                // outx negative -> flow from border to cell i,j-1
                newConcentrations(i, j - 1) -=
                    outx(i, j) * bc.getBoundaryElementValue(BC_SIDE_LEFT, i);
              }
            }
          }
          // flow between inner cells
          else {
            // outx positive -> flow from cell i,j-1 towards cell i,j
            if (outx(i, j) > 0) {
              newConcentrations(i, j - 1) -=
                  outx(i, j) * oldConcentrations(i, j - 1);
              newConcentrations(i, j) +=
                  outx(i, j) * oldConcentrations(i, j - 1);
            }
            // outx negative -> flow from cell i,j toward cell i,j-1
            else if (outx(i, j) < 0) {
              newConcentrations(i, j - 1) -=
                  outx(i, j) * oldConcentrations(i, j);
              newConcentrations(i, j) += outx(i, j) * oldConcentrations(i, j);
            }
          }
        }
      }
// calculate sum in y-direction
// parallelize outer loop over columns to ensure thread-safety, each thread only
// modifies cells within its column
#pragma omp parallel for num_threads(numThreads)
      for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows + 1; i++) {
          if (i == 0) {
            if (bc.getBoundaryElementType(BC_SIDE_TOP, j) != BC_TYPE_CLOSED) {
              if (outy(i, j) > 0) {
                // outy positive -> flow from border to cell i,j
                newConcentrations(i, j) +=
                    outy(i, j) * bc.getBoundaryElementValue(BC_SIDE_TOP, j);
              } else if (outy(i, j) < 0) {
                // outy negative -> flow from i,j towards border
                newConcentrations(i, j) += outy(i, j) * oldConcentrations(i, j);
              }
            }
          } else if (i == rows) {
            if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, j) !=
                BC_TYPE_CLOSED) {
              if (outy(i, j) > 0) {
                // outy positive-> flow from i-1,j towards border
                newConcentrations(i - 1, j) -=
                    outy(i, j) * oldConcentrations(i - 1, j);
              } else if (outy(i, j) < 0) {
                // outy negative -> flow from border to cell i,j-1
                newConcentrations(i - 1, j) -=
                    outy(i, j) * bc.getBoundaryElementValue(BC_SIDE_BOTTOM, j);
              }
            }
          }
          // flow between inner cells
          else {
            // outy positive -> flow from cell i-1,j towards cell i,j
            if (outy(i, j) > 0) {
              newConcentrations(i - 1, j) -=
                  outy(i, j) * oldConcentrations(i - 1, j);
              newConcentrations(i, j) +=
                  outy(i, j) * oldConcentrations(i - 1, j);
            }
            // outy negative -> flow from cell i,j toward cell i-1,j
            else if (outy(i, j) < 0) {
              newConcentrations(i - 1, j) -=
                  outy(i, j) * oldConcentrations(i, j);
              newConcentrations(i, j) += outy(i, j) * oldConcentrations(i, j);
            }
          }
        }
      }
      newConcentrations =
          oldConcentrations +
          newConcentrations * (innerTimestep / (porosity * volume));
    }
  }

  void run() {
    std::string filename;
    if (csv_output >= CSV_OUTPUT_ON) {
      filename = createCSVfile("Concentrations");
    }

    if (Steady == false) {
      velocities.setTimestep(timestep);
      velocities.setIterations(1);
    }

    for (int i = 0; i < iterations; i++) {
      if (csv_output >= CSV_OUTPUT_VERBOSE) {
        printConcentrationCSV(filename);
      }
      // if steady==false update charge and velocities with equal timestep
      if (Steady == false) {
        velocities.run();
      }
      adv();
    }

    if (csv_output >= CSV_OUTPUT_ON) {
      printConcentrationCSV(filename);
    }
  }

private:
  Grid<T> &grid;
  Boundary<T> &bc;
  Velocities<T> &velocities;
  bool Steady{true};
  int iterations{-1};
  int innerIterations{1};
  T timestep{-1};
  int numThreads{omp_get_num_procs()};
  T porosity{1};
  const Eigen::MatrixX<T> &outx;
  const Eigen::MatrixX<T> &outy;
  CSV_OUTPUT csv_output{CSV_OUTPUT_OFF};
};
} // namespace tug
