/**
 * @file Advection.hpp
 * @brief API of Advection class, holding information for a simulation of
 * advection in 1D and 2D grids.
 */

#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <Eigen/Core>
#include <tug/Advection/Velocities.hpp>
#include <tug/Boundary.hpp>
#include <tug/Core/BaseSimulation.hpp>
#include <tug/Core/Matrix.hpp>
#include <tug/Core/TugUtils.hpp>

namespace tug {

template <class T, HYDRAULIC_MODE hyd_mode = HYDRAULIC_MODE::STEADY_STATE,
          HYDRAULIC_RESOLVE hyd_resolve = HYDRAULIC_RESOLVE::IMPLICIT>
class Advection : public BaseSimulationGrid<T> {
private:
  T timestep{-1};
  int numThreads{omp_get_num_procs()};

  Velocities<T, hyd_mode, hyd_resolve> &velocities;
  RowMajMat<T> porosity;

  void init_advection() {
    tug_assert(this->rows() == velocities.rows() &&
                   this->cols() == velocities.cols(),
               "Advection grid and Velocities must have the same dimensions");
    this->porosity =
        RowMajMat<T>::Constant(this->rows(), this->cols(), static_cast<T>(1));
  }

public:
  /**
   * @brief Construct a new Advection object from an Eigen RowMajMat (1D or 2D).
   *
   * @param origin Concentration matrix.
   * @param _velocities Velocities simulation object of matching dimensions.
   */
  Advection(RowMajMat<T> &origin,
            Velocities<T, hyd_mode, hyd_resolve> &_velocities)
      : BaseSimulationGrid<T>(origin), velocities(_velocities) {
    init_advection();
  }

  /**
   * @brief Construct a new 2D Advection object from a raw data pointer.
   *
   * @param data Pointer to concentration array.
   * @param rows Number of rows.
   * @param cols Number of columns.
   * @param _velocities Velocities simulation object of matching dimensions.
   */
  Advection(T *data, std::size_t rows, std::size_t cols,
            Velocities<T, hyd_mode, hyd_resolve> &_velocities)
      : BaseSimulationGrid<T>(data, rows, cols), velocities(_velocities) {
    init_advection();
  }

  /**
   * @brief Construct a new 1D Advection object from a raw data pointer.
   *
   * @param data Pointer to concentration array.
   * @param length Length of the 1D grid.
   * @param _velocities Velocities simulation object of matching dimensions.
   */
  Advection(T *data, std::size_t length,
            Velocities<T, hyd_mode, hyd_resolve> &_velocities)
      : BaseSimulationGrid<T>(data, 1, length), velocities(_velocities) {
    init_advection();
  }

  /**
   * @brief Sets the porosity matrix of the medium.
   *
   * @param porosity New porosity matrix with values in (0, 1].
   */
  void setPorosity(const RowMajMat<T> &porosity) {
    tug_assert(porosity.rows() == this->rows() &&
                   porosity.cols() == this->cols(),
               "Porosity matrix must have the same dimensions as the grid");
    tug_assert(porosity.minCoeff() > 0 && porosity.maxCoeff() <= 1,
               "Porosity must be greater than 0 and at most 1");

    this->porosity = porosity;
  }

  /**
   * @brief Sets a uniform scalar porosity across the entire grid.
   *
   * @param constantPorosity Scalar porosity value in (0, 1].
   */
  void setPorosity(T constantPorosity) {
    tug_assert(constantPorosity > 0 && constantPorosity <= 1,
               "Porosity must be greater than 0 and at most 1");

    this->porosity =
        RowMajMat<T>::Constant(this->rows(), this->cols(), constantPorosity);
  }

  /**
   * @brief Returns the porosity matrix.
   *
   * @return const RowMajMat<T>& Reference to the porosity matrix.
   */
  const RowMajMat<T> &getPorosity() const { return this->porosity; }

  /**
   * @brief Set the size of the timestep. Must be greater than zero.
   *
   * @param timestep Size of the timestep.
   */
  void setTimestep(T timestep) override {
    tug_assert(timestep > 0, "Timestep must be greater than zero");
    this->timestep = timestep;
  }

  /**
   * @brief Returns the currently configured timestep.
   *
   * @return T Current timestep.
   */
  T getTimestep() const { return this->timestep; }

  /**
   * @brief Set the number of desired OpenMP threads.
   *
   * @param num_threads Number of desired threads (1 .. omp_get_num_procs()).
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
   * @brief Outputs the current concentrations of the grid on the console.
   */
  void printConcentrationsConsole() const {
    std::cout << this->getConcentrationMatrix() << std::endl;
    std::cout << std::endl;
  }

  /**
   * @brief Creates a CSV file with a name containing the current simulation
   * parameters.
   *
   * @return std::string Filename with configured simulation parameters.
   */
  std::string createCSVfile() const {
    std::ofstream file;
    int appendIdent = 0;

    std::string row = std::to_string(this->rows());
    std::string col = std::to_string(this->cols());
    std::string numIterations = std::to_string(this->getIterations());

    std::string filename =
        "ADV_" + row + "_" + col + "_" + numIterations + ".csv";

    while (std::filesystem::exists(filename)) {
      appendIdent += 1;
      filename = "ADV_" + row + "_" + col + "_" + numIterations + "-" +
                 std::to_string(appendIdent) + ".csv";
    }

    file.open(filename);
    if (!file) {
      throw std::runtime_error("Failed to open file: " + filename);
    }

    if (this->getOutputCSV() == CSV_OUTPUT::XTREME) {
      const auto &bc = this->getBoundaryConditions();
      Eigen::IOFormat one_row(Eigen::StreamPrecision, Eigen::DontAlignCols, "",
                              " ");
      file << bc.getBoundarySideValues(BC_SIDE_LEFT).format(one_row)
           << std::endl;
      file << bc.getBoundarySideValues(BC_SIDE_RIGHT).format(one_row)
           << std::endl;
      if (this->getDim() == 2) {
        file << bc.getBoundarySideValues(BC_SIDE_TOP).format(one_row)
             << std::endl;
        file << bc.getBoundarySideValues(BC_SIDE_BOTTOM).format(one_row)
             << std::endl;
      }
      file << std::endl << std::endl;
    }

    file.close();
    return filename;
  }

  /**
   * @brief Writes the currently calculated concentration values of the grid
   * into the CSV file with the passed filename.
   *
   * @param filename Name of the file to which the concentration values are
   * written.
   */
  void printConcentrationsCSV(const std::string &filename) const {
    std::ofstream file;
    file.open(filename, std::ios_base::app);
    if (!file) {
      throw std::runtime_error("Failed to open file: " + filename);
    }

    Eigen::IOFormat do_not_align(Eigen::StreamPrecision, Eigen::DontAlignCols);
    file << this->getConcentrationMatrix().format(do_not_align) << std::endl;
    file << std::endl << std::endl;
    file.close();
  }

  /**
   * @brief Starts the advection simulation with the configured parameters.
   */
  void run() override {
    tug_assert(this->getTimestep() > 0, "Timestep is not set!");
    tug_assert(this->getIterations() > 0, "Number of iterations are not set!");

    if (this->getDim() == 1) {
      this->setDomain(velocities.domainX());
    } else {
      this->setDomain(velocities.domainY(), velocities.domainX());
    }

    std::string filename;
    if (this->getOutputConsole() > CONSOLE_OUTPUT::OFF) {
      printConcentrationsConsole();
    }
    if (this->getOutputCSV() > CSV_OUTPUT::OFF) {
      filename = createCSVfile();
    }

    auto begin = std::chrono::high_resolution_clock::now();

    this->applyInnerBoundaries();

    if constexpr (hyd_mode == HYDRAULIC_MODE::STEADY_STATE) {
      velocities.run();
    } else {
      if (velocities.getTimestep() <= 0) {
        velocities.setTimestep(this->timestep);
      }
    }

    for (int i = 0; i < this->getIterations(); i++) {
      if (this->getOutputConsole() == CONSOLE_OUTPUT::VERBOSE && i > 0) {
        printConcentrationsConsole();
      }
      if (this->getOutputCSV() >= CSV_OUTPUT::VERBOSE) {
        printConcentrationsCSV(filename);
      }

      if constexpr (hyd_mode == HYDRAULIC_MODE::TRANSIENT) {
        velocities.run();
      }
      startAdvection();
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
      const std::string dimString = std::to_string(this->getDim()) + "D";
      std::cout << "ADV" << dimString << ":: run() finished in "
                << milliseconds.count() << "ms" << std::endl;
    }
  }

private:
  /**
   * @brief Calculates material transport for one timestep using a finite-volume
   * first-order upwind scheme with automatic CFL substepping.
   */
  void startAdvection() {
    const std::size_t rows = this->rows();
    const std::size_t cols = this->cols();
    const T dx = this->deltaCol();
    const T dy = (this->getDim() == 2) ? this->deltaRow() : static_cast<T>(1);

    RowMajMatMap<T> &newConcentrations = this->getConcentrationMatrix();
    const RowMajMat<T> &veloX = this->velocities.getVelocitiesX();
    const Boundary<T> &bc = this->getBoundaryConditions();

    // Calculate Courant-Friedrichs-Lewy condition
    const T maxVx = veloX.cwiseAbs().maxCoeff() / dx;
    T maxRate = maxVx;
    if (this->getDim() == 2) {
      const RowMajMat<T> &veloY = this->velocities.getVelocitiesY();
      const T maxVy = veloY.cwiseAbs().maxCoeff() / dy;
      maxRate += maxVy;
    }

    if (maxRate == 0) {
      return;
    }

    const T cfl = this->porosity.minCoeff() / maxRate;
    const int innerSteps =
        static_cast<int>(std::ceil(this->timestep / cfl));
    const T innerTimestep = this->timestep / static_cast<T>(innerSteps);

    if (this->getOutputConsole() > CONSOLE_OUTPUT::OFF) {
      std::cout << "CFL (adv) timestep: " << cfl << std::endl;
      std::cout << "Inner iterations (adv): " << innerSteps << std::endl;
      std::cout << "Inner timestep (adv): " << innerTimestep << std::endl;
    }

    RowMajMat<T> fluxX(rows, cols + 1);
    RowMajMat<T> fluxY;
    if (this->getDim() == 2) {
      fluxY.resize(rows + 1, cols);
    }

    for (int k = 0; k < innerSteps; k++) {
      const RowMajMat<T> oldConcentrations = newConcentrations;

      // Update flux with concentrations in x-direction
#pragma omp parallel for num_threads(numThreads)
      for (std::size_t row_i = 0; row_i < rows; row_i++) {
        for (std::size_t col_i = 0; col_i <= cols; col_i++) {
          T &currentFlux = fluxX(row_i, col_i);
          const T currentVelo = veloX(row_i, col_i);

          if (col_i == 0) {
            const auto bcLeft = bc.getBoundaryElement(BC_SIDE_LEFT, row_i);
            if (bcLeft.getType() == BC_TYPE_CLOSED) {
              currentFlux = 0;
            } else if (currentVelo > 0) {
              currentFlux = currentVelo * bcLeft.getValue();
            } else {
              currentFlux = currentVelo * oldConcentrations(row_i, 0);
            }
          } else if (col_i == cols) {
            const auto bcRight = bc.getBoundaryElement(BC_SIDE_RIGHT, row_i);
            if (bcRight.getType() == BC_TYPE_CLOSED) {
              currentFlux = 0;
            } else if (currentVelo > 0) {
              currentFlux = currentVelo * oldConcentrations(row_i, cols - 1);
            } else {
              currentFlux = currentVelo * bcRight.getValue();
            }
          } else {
            const std::size_t cellIndex =
                (currentVelo > 0) ? (col_i - 1) : col_i;
            currentFlux = currentVelo * oldConcentrations(row_i, cellIndex);
          }
        }
      }

      // Update flux with concentrations in y-direction (2D only)
      if (this->getDim() == 2) {
        const RowMajMat<T> &veloY = this->velocities.getVelocitiesY();
#pragma omp parallel for num_threads(numThreads)
        for (std::size_t row_i = 0; row_i <= rows; row_i++) {
          for (std::size_t col_i = 0; col_i < cols; col_i++) {
            T &currentFlux = fluxY(row_i, col_i);
            const T currentVelo = veloY(row_i, col_i);

            if (row_i == 0) {
              const auto bcTop = bc.getBoundaryElement(BC_SIDE_TOP, col_i);
              if (bcTop.getType() == BC_TYPE_CLOSED) {
                currentFlux = 0;
              } else if (currentVelo > 0) {
                currentFlux = currentVelo * bcTop.getValue();
              } else {
                currentFlux = currentVelo * oldConcentrations(0, col_i);
              }
            } else if (row_i == rows) {
              const auto bcBottom =
                  bc.getBoundaryElement(BC_SIDE_BOTTOM, col_i);
              if (bcBottom.getType() == BC_TYPE_CLOSED) {
                currentFlux = 0;
              } else if (currentVelo > 0) {
                currentFlux = currentVelo * oldConcentrations(rows - 1, col_i);
              } else {
                currentFlux = currentVelo * bcBottom.getValue();
              }
            } else {
              const std::size_t cellIndex =
                  (currentVelo > 0) ? (row_i - 1) : row_i;
              currentFlux = currentVelo * oldConcentrations(cellIndex, col_i);
            }
          }
        }
      }

      // Update concentrations
#pragma omp parallel for num_threads(numThreads)
      for (std::size_t row_i = 0; row_i < rows; row_i++) {
        for (std::size_t col_i = 0; col_i < cols; col_i++) {
          const bool isInnerBC =
              (this->getDim() == 1)
                  ? bc.getInnerBoundary(col_i).first
                  : bc.getInnerBoundary(row_i, col_i).first;
          if (isInnerBC) {
            continue;
          }

          T totalRate = (fluxX(row_i, col_i) - fluxX(row_i, col_i + 1)) / dx;
          if (this->getDim() == 2) {
            totalRate += (fluxY(row_i, col_i) - fluxY(row_i + 1, col_i)) / dy;
          }

          newConcentrations(row_i, col_i) =
              oldConcentrations(row_i, col_i) +
              (innerTimestep / this->porosity(row_i, col_i)) * totalRate;
        }
      }
    }
  }
};
} // namespace tug
