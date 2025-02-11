/**
 * @file Advection.hpp
 * @brief API of Advection class, holding information for a simulation of
 * advection. Holds a predifined Grid object, Boundary object and Velocities
 * object
 */

#pragma once

#include "tug/Core/Matrix.hpp"
#include <Eigen/src/Core/Array.h>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <tug/Boundary.hpp>

#include <Eigen/Eigen>
#include <tug/Core/Numeric/BTCS.hpp>
#include <tug/Core/Numeric/FTCS.hpp>
#include <tug/Core/TugUtils.hpp>
#include <tug/Diffusion/Diffusion.hpp>

#include <tug/Advection/Velocities.hpp>

using namespace Eigen;
namespace tug {
template <class T, HYDRAULIC_MODE hyd_mode, HYDRAULIC_RESOLVE hyd_resolve>
class Advection : public BaseSimulationGrid<T> {
private:
  T timestep{-1};
  int numThreads{omp_get_num_procs()};

  Velocities<T, hyd_mode, hyd_resolve> &velocities;
  RowMajMat<T> porosity;

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
  Advection(RowMajMat<T> &origin,
            Velocities<T, hyd_mode, hyd_resolve> &_velocities)
      : BaseSimulationGrid<T>(origin), velocities(_velocities) {
    tug_assert(origin.rows() == velocities.getConcentrationMatrix().rows() &&
                   origin.cols() == velocities.getConcentrationMatrix().cols(),
               "Advection grid and Velocities must have the same dimensions");
  };

  Advection(T *data, std::size_t rows, std::size_t cols,
            Velocities<T, hyd_mode, hyd_resolve> &_velocities)
      : BaseSimulationGrid<T>(data, rows, cols), velocities(_velocities) {
    tug_assert(rows == velocities.getConcentrationMatrix().rows() &&
                   cols == velocities.getConcentrationMatrix().cols(),
               "Advection grid and Velocities must have the same dimensions");
  };

  /**
   * @brief Sets the porosity of the medium
   *
   * @param porosity new porosity value
   */
  void setPorosity(const RowMajMat<T> &porosity) {
    tug_assert(porosity.rows() == this->rows() &&
                   porosity.cols() == this->cols(),
               "Porosity matrix must have the same dimensions as the grid");
    tug_assert(porosity.minCoeff() >= 0 && porosity.maxCoeff() <= 1,
               "Porosity must be a value between 0 and 1 (inclusive)");

    this->porosity = porosity;
  }

  /**
   * @brief Set the size of the timestep. Must be greater than zero
   *
   * @param timestep Size of the timestep
   */
  void setTimestep(T timestep) {
    tug_assert(timestep > 0, "Timestep must be greater than zero");
    this->timestep = timestep;
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

  void run() {
    this->setDomain(velocities.domainX(), velocities.domainY());

    this->applyInnerBoundaries();

    if constexpr (hyd_mode == HYDRAULIC_MODE::STEADY_STATE) {
      velocities.run();
    }

    for (int i = 0; i < this->getIterations(); i++) {
      if constexpr (hyd_mode == HYDRAULIC_MODE::TRANSIENT) {
        velocities.run();
      }
      startAdvection();
    }
  }

private:
  /**
   * @brief function calculating material transport for one timestep
   */
  void startAdvection() {
    const std::size_t rows = this->rows();
    const std::size_t cols = this->cols();
    const T dx = this->deltaCol();
    const T dy = this->deltaRow();
    const T volume = this->deltaCol() * this->deltaRow();

    RowMajMatMap<T> &newConcentrations = this->getConcentrationMatrix();

    RowMajMat<T> veloX = this->velocities.getVelocitiesX();
    RowMajMat<T> veloY = this->velocities.getVelocitiesY();
    const Boundary<T> &bc = this->getBoundaryConditions();

    // Calculate Courant-Levy-Frederich condition
    const T maxFx = std::max(abs(veloX.maxCoeff()), abs(veloX.minCoeff()));
    const T maxFy = std::max(abs(veloY.maxCoeff()), abs(veloY.minCoeff()));
    const T maxFlux = std::max(maxFx, maxFy);

    tug_assert(maxFlux != 0, "Division by zero: maxF is zero.");

    // TODO: sustitute 1 by porisity
    const RowMajMat<T> porevolume =
        RowMajMat<T>::Constant(rows, cols, volume) * 1;

    const T cfl = (porevolume / maxFlux).minCoeff();
    const int innerSteps = (int)ceil(timestep / cfl);
    const T innerTimestep = timestep / innerSteps;

    std::cout << "CFL (adv) timestep: " << cfl << std::endl;
    std::cout << "Inner iterations (adv): " << innerSteps << std::endl;
    std::cout << "Inner timestep (adv): " << innerTimestep << std::endl;

    // const RowMajMat<T> multiplier = porevolume * (1 / innerTimestep);

    RowMajMat<T> fluxX(rows, cols + 1);
    RowMajMat<T> fluxY(rows + 1, cols);

    for (std::size_t k = 0; k < innerSteps; k++) {
      const RowMajMat<T> oldConcentrations = newConcentrations;

      // Update flux with concentrations in x-direction
      for (std::size_t row_i = 0; row_i < fluxX.rows(); row_i++) {
        for (std::size_t col_i = 0; col_i < fluxX.cols(); col_i++) {
          T &currentFlux = fluxX(row_i, col_i);
          const T &currentVelo = veloX(row_i, col_i);
          const std::int32_t cellIndex = (currentVelo > 0) ? col_i - 1 : col_i;

          if (cellIndex < 0 || cellIndex >= cols) {
            const auto bcElement = bc.getBoundaryElement(
                cellIndex < 0 ? BC_SIDE_LEFT : BC_SIDE_RIGHT, row_i);
            // if we have a constant boundary, we use the value of the boundary
            if (bcElement.getType() == BC_TYPE_CONSTANT) {
              currentFlux =
                  currentVelo * bcElement.getValue() * (innerTimestep) / dx;
              continue;
            }

            // otherwise it is a closed boundary
            currentFlux = 0;
            continue;
          }
          // otherwise we use the concentration of the inflow/outflow cell,
          // multiplied by the velocity
          currentFlux = currentVelo * oldConcentrations(row_i, cellIndex) *
                        innerTimestep / dx;
        }
      }

      // Update flux with concentrations in y-direction
      for (std::size_t row_i = 0; row_i < fluxY.rows(); row_i++) {
        for (std::size_t col_i = 0; col_i < fluxY.cols(); col_i++) {
          T &currentFlux = fluxY(row_i, col_i);
          const T &currentVelo = veloY(row_i, col_i);
          const std::int32_t cellIndex = (currentVelo > 0) ? row_i - 1 : row_i;

          if (cellIndex < 0 || cellIndex >= rows) {
            const auto bcElement = bc.getBoundaryElement(
                cellIndex < 0 ? BC_SIDE_TOP : BC_SIDE_BOTTOM, col_i);
            // if we have a constant boundary, we use the value of the
            // boundary
            if (bcElement.getType() == BC_TYPE_CONSTANT) {
              currentFlux =
                  currentVelo * bcElement.getValue() * (2 * innerTimestep) / dy;
              continue;
            }

            // otherwise it is a closed boundary
            currentFlux = 0;
            continue;
          }
          // otherwise we use the concentration of the inflow/outflow cell,
          // multiplied by the velocity
          currentFlux = currentVelo * oldConcentrations(cellIndex, col_i) *
                        innerTimestep / dy;
        }
      }

      // Update concentrations
      for (std::size_t row_i = 0; row_i < rows; row_i++) {
        for (std::size_t col_i = 0; col_i < cols; col_i++) {
          const T horizontalFlux =
              fluxX(row_i, col_i) - fluxX(row_i, col_i + 1);
          const T verticalFlux = fluxY(row_i, col_i) - fluxY(row_i + 1, col_i);
          newConcentrations(row_i, col_i) =
              oldConcentrations(row_i, col_i) + horizontalFlux + verticalFlux;
        }
      }

      // const Boundary<T> leftCorner = bc.getBoundaryElement(BC_SIDE_LEFT, 0);
      // newConcentrations(0, 0) = oldConcentrations(0, 0) + ;
      //
      // for (std::size_t col_i = 0; col_i < cols; col_i++) {
      //   T horizontalFlux = oldConcentrations(0, col_i);
      //   if (col_i == 0) {
      //     switch (bc.getBoundaryElementType(BC_SIDE_LEFT, 0)) {
      //     case BC_TYPE_CLOSED: {
      //       break;
      //     }
      //     case BC_TYPE_CONSTANT: {
      //       newConcentrations(0, 0) =
      //           oldConcentrations(0, 0) +
      //           (bc.getBoundaryElementValue(BC_SIDE_LEFT, 0) * veloX(0, 0)) *
      //               (2 * innerTimestep) / dx;
      //       newConcentrations(rows - 1, 0) =
      //           oldConcentrations(rows - 1, 0) +
      //           (bc.getBoundaryElementValue(BC_SIDE_LEFT, 0) *
      //            veloX(rows - 1, 0)) *
      //               (2 * innerTimestep) / dx;
      //       break;
      //     }
      //     }
      //     newConcentrations(0, 0) -=
      //         veloX(0, 1) * oldConcentrations(0, 1) * innerTimestep / dx;
      //   }
      // }
      //
      // for (std::size_t row_i = 1; row_i < rows - 1; row_i++) {
      //   for (std::size_t col_i = 1; col_i < cols - 1; col_i++) {
      //     const T fluxHorizontal =
      //         veloX(row_i, col_i) * oldConcentrations(row_i, col_i - 1) -
      //         veloX(row_i, col_i + 1) * oldConcentrations(row_i, col_i);
      //     const T fluxVertical =
      //         veloY(row_i, col_i) * oldConcentrations(row_i - 1, col_i) -
      //         veloY(row_i + 1, col_i) * oldConcentrations(row_i, col_i);
      //     newConcentrations(row_i, col_i) =
      //         oldConcentrations(row_i, col_i) +
      //         (fluxHorizontal / dx + fluxVertical / dy) * innerTimestep;
      //   }
      // }

      // Calculate sum of incoming/outgoing Flow*Concentration in x-direction in
      // each cell #pragma omp parallel for num_threads(numThreads)
      //       for (int i = 0; i < rows; i++) {
      //         for (int j = 0; j < cols + 1; j++) {
      //           if (j == 0) {
      //             if (bc.getBoundaryElementType(BC_SIDE_LEFT, i) !=
      //             BC_TYPE_CLOSED) {
      //               if (veloX(i, j) > 0) {
      //                 // outx positive -> flow from border to cell i,j
      //                 newConcentrations(i, j) +=
      //                     veloX(i, j) *
      //                     bc.getBoundaryElementValue(BC_SIDE_LEFT, i);
      //               } else if (veloX(i, j) < 0) {
      //                 // outx negative -> flow from i,j towards border
      //                 newConcentrations(i, j) += veloX(i, j) *
      //                 oldConcentrations(i, j);
      //               }
      //             }
      //           } else if (j == cols) {
      //             if (bc.getBoundaryElementType(BC_SIDE_RIGHT, i) !=
      //             BC_TYPE_CLOSED) {
      //               if (veloX(i, j) > 0) {
      //                 // outx positive-> flow from i,j-1 towards border
      //                 newConcentrations(i, j - 1) -=
      //                     veloX(i, j) * oldConcentrations(i, j - 1);
      //               } else if (veloX(i, j) < 0) {
      //                 // outx negative -> flow from border to cell i,j-1
      //                 newConcentrations(i, j - 1) -=
      //                     veloX(i, j) *
      //                     bc.getBoundaryElementValue(BC_SIDE_LEFT, i);
      //               }
      //             }
      //           }
      //           // flow between inner cells
      //           else {
      //             // outx positive -> flow from cell i,j-1 towards cell i,j
      //             if (veloX(i, j) > 0) {
      //               newConcentrations(i, j - 1) -=
      //                   veloX(i, j) * oldConcentrations(i, j - 1);
      //               newConcentrations(i, j) +=
      //                   veloX(i, j) * oldConcentrations(i, j - 1);
      //             }
      //             // outx negative -> flow from cell i,j toward cell i,j-1
      //             else if (veloX(i, j) < 0) {
      //               newConcentrations(i, j - 1) +=
      //                   veloX(i, j) * oldConcentrations(i, j);
      //               newConcentrations(i, j) += veloX(i, j) *
      //               oldConcentrations(i, j);
      //             }
      //           }
      //         }
      //       }
      // // calculate sum in y-direction
      // // parallelize outer loop over columns to ensure thread-safety, each
      // thread only
      // // modifies cells within its column
      // #pragma omp parallel for num_threads(numThreads)
      //       for (int j = 0; j < cols; j++) {
      //         for (int i = 0; i < rows + 1; i++) {
      //           if (i == 0) {
      //             if (bc.getBoundaryElementType(BC_SIDE_TOP, j) !=
      //             BC_TYPE_CLOSED) {
      //               if (veloY(i, j) > 0) {
      //                 // outy positive -> flow from border to cell i,j
      //                 newConcentrations(i, j) +=
      //                     veloY(i, j) *
      //                     bc.getBoundaryElementValue(BC_SIDE_TOP, j);
      //               } else if (veloY(i, j) < 0) {
      //                 // outy negative -> flow from i,j towards border
      //                 newConcentrations(i, j) += veloY(i, j) *
      //                 oldConcentrations(i, j);
      //               }
      //             }
      //           } else if (i == rows) {
      //             if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, j) !=
      //                 BC_TYPE_CLOSED) {
      //               if (veloY(i, j) > 0) {
      //                 // outy positive-> flow from i-1,j towards border
      //                 newConcentrations(i - 1, j) -=
      //                     veloY(i, j) * oldConcentrations(i - 1, j);
      //               } else if (veloY(i, j) < 0) {
      //                 // outy negative -> flow from border to cell i,j-1
      //                 newConcentrations(i - 1, j) -=
      //                     veloY(i, j) *
      //                     bc.getBoundaryElementValue(BC_SIDE_BOTTOM, j);
      //               }
      //             }
      //           }
      //           // flow between inner cells
      //           else {
      //             // outy positive -> flow from cell i-1,j towards cell i,j
      //             if (veloY(i, j) > 0) {
      //               newConcentrations(i - 1, j) -=
      //                   veloY(i, j) * oldConcentrations(i - 1, j);
      //               newConcentrations(i, j) +=
      //                   veloY(i, j) * oldConcentrations(i - 1, j);
      //             }
      //             // outy negative -> flow from cell i,j toward cell i-1,j
      //             else if (veloY(i, j) < 0) {
      //               newConcentrations(i - 1, j) -=
      //                   veloY(i, j) * oldConcentrations(i, j);
      //               newConcentrations(i, j) += veloY(i, j) *
      //               oldConcentrations(i, j);
      //             }
      //           }
      //         }
      //       }

      // for (std::size_t row_i = 0; row_i < rows; row_i++) {
      //   for (std::size_t col_i = 0; col_i < cols; col_i++) {
      //     newConcentrations(row_i, col_i) =
      //         oldConcentrations(row_i, col_i) +
      //         newConcentrations(row_i, col_i) * multiplier(row_i, col_i);
      //   }
      // }
    }
  }
};
} // namespace tug
