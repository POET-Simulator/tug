/**
 * @file Velocities.hpp
 * @brief API of Velocities class, holding information for a simulation of
 * Hydraulic Charge and Darcy-Velocities. Holds a predifined Grid object and
 * Boundary object.
 *
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
#include <vector>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <tug/Core/Numeric/BTCS.hpp>
#include <tug/Core/Numeric/FTCS.hpp>
#include <tug/Core/TugUtils.hpp>
#include <tug/Diffusion.hpp>

using namespace Eigen;
namespace tug {
template <class T> class Velocities {
public:
  /**
   * @brief Construct a new Velocities object, used to calculate Hydraulic
   * Charge and Darcy-Velocities. A timestep and a number of iterations can be
   * set. By default iterations is set to -1. If the number of iterations is set
   * to a value below 1 the simulation will run until the Hydraulic Charge
   * converges. The Epsilon value to check convergence can be set, the default
   * is 1E-5. CSV Output is off by default.
   *
   * @param grid Valid grid object
   * @param bc Valid boundary condition object
   */
  Velocities(Grid<T> &_grid, Boundary<T> &_bc) : grid(_grid), bc(_bc) {
    outx = MatrixX<T>::Constant(grid.getRow(), grid.getCol() + 1, 0);
    outy = MatrixX<T>::Constant(grid.getRow() + 1, grid.getCol(), 0);
    center = std::make_pair(grid.getRow() / 2, grid.getCol() / 2);
  };

  /**
   * @brief Sets a fixed, constant hydraulic charge at domain center.
   *
   * @param inj_h fixed hydraulic charge at domain center.
   */
  void setInjh(T inj_h) {
    if (inj_h >= 0) {
      this->inj_h = inj_h;
      RowMajMat<T> &concentrations = grid.getConcentrations();
      concentrations(center.first, center.second) = inj_h;
      injhIsSet = true;
    } else {
      throw std::invalid_argument("Fixed hydraulic charge can not be negative");
    }
  };

  /**
   * @brief Sets a constant permeability coefficient
   * @param K constant permeability coefficient
   */
  void setK(T K) {
    this->K = K;
    MatrixXd alphax = MatrixXd::Constant(grid.getRow(), grid.getCol(), K);
    MatrixXd alphay = MatrixXd::Constant(grid.getRow(), grid.getCol(), K);
    grid.setAlpha(alphax, alphay);
  };

  /**
   * @brief Set the epsilon value, the relativ error allowed for convergence
   *
   * @param epsilon the new epsilon value
   */
  void setEpsilon(T epsilon) {
    if (0 <= epsilon && epsilon < 1) {
      this->epsilon = epsilon;
    } else {
      throw std::invalid_argument(
          "Relative Error epsilon must be between 0 and 1");
    }
  }

  /**
   * @brief Set the timestep per iteration
   *
   * @param timestep timestep per iteration
   */
  void setTimestep(T timestep) {
    if (timestep <= 0) {
      throw std::invalid_argument("Timestep must be greater than zero");
    }
    this->timestep = timestep;
    const T deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();
    const T deltaRowSquare = grid.getDeltaRow() * grid.getDeltaRow();
    const T minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);
    T cfl = minDeltaSquare / (4 * K);
    if (timestep > cfl) {
      this->innerIterations = (int)ceil(timestep / cfl);
      this->timestep = timestep / (double)innerIterations;
      std::cerr << "Warning :: Timestep was adjusted, because of stability "
                   "conditions. Time duration was approximately preserved by "
                   "adjusting internal number of iterations."
                << std::endl;
      std::cout << "FTCS" << "_" << "2D" << " :: Required "
                << this->innerIterations
                << " inner iterations with dt=" << this->timestep << std::endl;
    } else {
      this->innerIterations = 1;
    }
  };

  /**
   * @brief Set the number of iterations. If set to a number smaller than 1,
   * calculation will terminate at convergence
   *
   * @param iterations Number of desired iterations
   */

  void setIterations(int iterations) { this->iterations = iterations; }

  /**
   * @brief Set the number of desired openMP Threads.
   *
   * @param num_threads Number of desired threads. Must have a value between
   *                    1 and the maximum available number of processors. The
   * maximum number of processors is set as the default case during Velocities
   * construction.
   */
  void setNumberThreads(int num_threads) {
    if (num_threads > 0 && num_threads <= omp_get_num_procs()) {
      this->numThreads = num_threads;
    } else {
      int maxThreadNumber = omp_get_num_procs();
      throw std::invalid_argument(
          "Number of threads exceeds the number of processor cores (" +
          std::to_string(maxThreadNumber) + ") or is less than 1.");
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
   *                       charge and velocities matrizes
   *                     - CSV_OUTPUT_VERBOSE: produce csv output with all
   *                       charge matrizes and and last velocities matrix
   */
  void setOutputCSV(CSV_OUTPUT csv_output) {
    if (csv_output < CSV_OUTPUT_OFF && csv_output > CSV_OUTPUT_VERBOSE) {
      throw std::invalid_argument("Invalid CSV output option given!");
    }
    this->csv_output = csv_output;
    if (csv_output >= CSV_OUTPUT_ON) {
      filename1 = createCSVfile("Charge");
      filename2 = createCSVfile("outx");
      filename3 = createCSVfile("outy");
    }
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
      exit(1);
    }

    file.close();

    return filename;
  }

  /**
   * @brief Writes the currently calculated Hydraulic Charge values of the grid
   *        into the CSV file with the passed filename.
   *
   * @param filename Name of the file to which the Hydraulic Charge values are
   *                 to be written.
   */
  void printChargeCSV(const std::string &filename) const {
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
   * @brief Reads a matrix stored in a CSV file and uses it for values of
   * Hydraulic Heads, Matrix and grid must be of equal size
   *
   * @param filename name of the CSV file
   */
  void readChargeCSV(std::string filename) {
    std::ifstream file(filename);
    std::string line;
    Eigen::MatrixXd matrix;

    if (file.is_open()) {
      while (std::getline(file, line)) {
        std::istringstream iss(line);
        double value;
        std::vector<double> row;
        while (iss >> value) {
          row.push_back(value);
        }
        if (!row.empty()) {
          if (matrix.rows() == 0) {
            matrix.resize(1, row.size());
          } else {
            matrix.conservativeResize(matrix.rows() + 1, Eigen::NoChange);
          }
          matrix.row(matrix.rows() - 1) =
              Eigen::VectorXd::Map(row.data(), row.size());
        }
      }
      file.close();
    } else {
      std::cerr << "Unable to open file: " << filename << std::endl;
    }
    if (matrix.rows() == grid.getRow() && matrix.cols() == grid.getCol()) {
      grid.setConcentrations(matrix);
      velocities();
      if (csv_output > CSV_OUTPUT_OFF) {
        printChargeCSV(filename1);
        printVelocitiesCSV(filename2, filename3);
      }
    } else {
      std::cerr << "gridsize and size of stored matrix dont align\n";
    }
  }

  /**
   * @brief Writes the current Darcy-velocities into a CSV file
   *
   * @param filenamex Name of the file to which velocities in direction x are
   * written
   * @param filenamey Name of the file to which velocities in direction y are
   * written
   */
  void printVelocitiesCSV(const std::string &filenamex,
                          const std::string &filenamey) const {
    std::ofstream filex;
    std::ofstream filey;

    filex.open(filenamex, std::ios_base::app);
    if (!filex) {
      exit(1);
    }
    filey.open(filenamey, std::ios_base::app);
    if (!filey) {
      exit(1);
    }

    Eigen::IOFormat do_not_align(Eigen::StreamPrecision, Eigen::DontAlignCols);
    filex << outx.format(do_not_align) << std::endl;
    filex << std::endl << std::endl;
    filex.close();

    filey << outy.format(do_not_align) << std::endl;
    filey << std::endl << std::endl;
    filey.close();
  }

  /**
   * @brief Calculate the new hydraulic charge using FTCS
   */
  void hydraulic_charge() {
    FTCS_2D(this->grid, this->bc, this->timestep, this->numThreads);
    if (injhIsSet == true) {
      RowMajMat<T> &concentrations = grid.getConcentrations();
      concentrations(center.first, center.second) = inj_h;
    }
  };

  /**
   * @brief checks if the matrix of Hydraulic Heads has converged
   *
   * @return bool true if for all corresponding cells of the matrices,
   * containing old and new Charge values, the relative error is below the
   * selected Epsilon
   */
  bool checkConvergance(Eigen::MatrixX<T> oldHeads,
                        Eigen::MatrixX<T> newHeads) {
    for (int i = 0; i < grid.getRow(); i++) {
      for (int j = 0; j < grid.getCol(); j++) {
        if (newHeads(i, j) != 0) {
          if (abs((oldHeads(i, j) - newHeads(i, j)) / newHeads(i, j)) >
              epsilon) {
            return false;
          }
        }
      }
    }
    return true;
  }

  /**
   * @brief Update the matrices containing Darcy velocities in x and y
   * directions
   */
  void velocities() {
    int rows = grid.getRow();
    int cols = grid.getCol();
    float dx = grid.getDeltaRow();
    float dy = grid.getDeltaCol();
    Eigen::MatrixX<T> concentrations = grid.getConcentrations();
// calculate outx
#pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols + 1; j++) {
        if (j == 0) {
          if (bc.getBoundaryElementType(BC_SIDE_LEFT, i) == BC_TYPE_CLOSED) {
            outx(i, j) = 0;
          } else {
            outx(i, j) = -K *
                         (concentrations(i, j) -
                          bc.getBoundaryElementValue(BC_SIDE_LEFT, i)) /
                         (dx / 2);
          }
        } else if (j == cols) {
          if (bc.getBoundaryElementType(BC_SIDE_RIGHT, i) == BC_TYPE_CLOSED) {
            outx(i, j) = 0;
          } else {
            outx(i, j) = -K *
                         (bc.getBoundaryElementValue(BC_SIDE_RIGHT, i) -
                          concentrations(i, j - 1)) /
                         (dx / 2);
          }

        } else {
          outx(i, j) =
              -K * (concentrations(i, j) - concentrations(i, j - 1)) / dx;
        }
      }
    }
// calculate outy
#pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < rows + 1; i++) {
      for (int j = 0; j < cols; j++) {
        if (i == 0) {
          if (bc.getBoundaryElementType(BC_SIDE_TOP, j) == BC_TYPE_CLOSED) {
            outy(i, j) = 0;
          } else {
            outy(i, j) = -K *
                         (concentrations(i, j) -
                          bc.getBoundaryElementValue(BC_SIDE_TOP, j)) /
                         (dy / 2);
          }
        } else if (i == rows) {
          if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, j) == BC_TYPE_CLOSED) {
            outy(i, j) = 0;
          } else {
            outy(i, j) = -K *
                         (bc.getBoundaryElementValue(BC_SIDE_BOTTOM, j) -
                          concentrations(i - 1, j)) /
                         (dy / 2);
          }
        } else {
          outy(i, j) =
              -K * (concentrations(i, j) - concentrations(i - 1, j)) / dy;
        }
      }
    }
  };

  /**
   * @brief Getter function for outx, the matrix containing velocities in
   * x-Direction; returns a reference to outx
   *
   * */
  const Eigen::MatrixX<T> &getOutx() const { return outx; }

  /**
   * @brief Getter function for outy, the matrix containing velocities in
   * y-Direction; return a reference to outy
   */
  const Eigen::MatrixX<T> &getOuty() const { return outy; }

  /**
   * @brief Simulation of hydraulic charge either until convergence,
   * or for a number of selected timesteps. Calculation of Darcy-velocities.
   */
  void run() {
    // if iterations < 1 calculate hydraulic charge until steady state is
    // reached
    if (iterations < 1) {
      // Calculate largest possible timestep, depending on K and gridsize
      const T deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();
      const T deltaRowSquare = grid.getDeltaRow() * grid.getDeltaRow();
      const T minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);
      setTimestep(minDeltaSquare / (4 * K));

      Eigen::MatrixX<T> oldConcentrations;
      do {
        oldConcentrations = grid.getConcentrations();
        if (csv_output >= CSV_OUTPUT_VERBOSE) {
          printChargeCSV(filename1);
        }
        for (int i = 0; i < (grid.getRow() + grid.getCol() - 2); i++) {
          hydraulic_charge();
        }
      } while (checkConvergance(oldConcentrations, grid.getConcentrations()) ==
               false);
    }
    // if iterations >= 1 calculate hydraulice charge for a given number of
    // iterations
    else {
      if (timestep == -1) {
        throw_invalid_argument("Timestep is not set");
      }
      for (int i = 0; i < iterations * innerIterations; i++) {
        hydraulic_charge();
      }
    }

    velocities();

    if (csv_output > CSV_OUTPUT_OFF) {
      printChargeCSV(filename1);
      printVelocitiesCSV(filename2, filename3);
    }
  };

private:
  int iterations{-1};
  int innerIterations{1};
  bool injhIsSet{false};
  T timestep{-1};
  T inj_h{1};
  T K{1};
  T epsilon{1E-5};
  int numThreads{omp_get_num_procs()};
  Grid<T> &grid;
  Boundary<T> &bc;
  CSV_OUTPUT csv_output{CSV_OUTPUT_OFF};
  Eigen::MatrixX<T> outx;
  Eigen::MatrixX<T> outy;
  std::pair<int, int> center;
  std::string filename1;
  std::string filename2;
  std::string filename3;
};
} // namespace tug
