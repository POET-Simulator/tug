/**
 * @file Velocities.hpp
 * @brief API of Velocities class, holding information for a simulation of
 * Hydraulic Charge and Darcy-Velocities. Holds a predifined Grid object and
 * Boundary object.
 *
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <string>

#include <tug/Boundary.hpp>
#include <tug/Core/Matrix.hpp>
#include <tug/Core/Numeric/BTCS.hpp>
#include <tug/Core/Numeric/FTCS.hpp>
#include <tug/Core/TugUtils.hpp>
#include <tug/Diffusion/Diffusion.hpp>
#include <tug/Grid.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs() 1
#endif

using namespace Eigen;
namespace tug {

enum HYDRAULIC_MODE { TRANSIENT, STEADY_STATE };
enum HYDRAULIC_RESOLVE { EXPLICIT, IMPLICIT };

template <class T, HYDRAULIC_MODE hyd_mode, HYDRAULIC_RESOLVE hyd_resolve>
class Velocities : public BaseSimulation {
private:
  int innerIterations{1};
  T timestep{-1};
  T epsilon{1E-5};
  int numThreads{omp_get_num_procs()};

  Grid<T> &grid;
  Boundary<T> &bc;

  RowMajMat<T> velocitiesX;
  RowMajMat<T> velocitiesY;

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
  Velocities(Grid<T> &_grid, Boundary<T> &_bc)
      : grid(_grid), bc(_bc), velocitiesX(grid.getRow(), grid.getCol() + 1),
        velocitiesY(grid.getRow() + 1, grid.getCol()) {};

  /**
   * @brief Set the epsilon value, the relativ error allowed for convergence
   *
   * @param epsilon the new epsilon value
   */
  void setEpsilon(T epsilon) {
    tug_assert(0 <= epsilon && epsilon < 1,
               "Relative Error epsilon must be between 0 and 1");

    this->epsilon = epsilon;
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

    const T maxK =
        std::max(grid.getAlphaX().maxCoeff(), grid.getAlphaY().maxCoeff());

    T cfl = minDeltaSquare / (4 * maxK);
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
   * @brief Getter function for outx, the matrix containing velocities in
   * x-Direction; returns a reference to outx
   *
   * */
  const RowMajMat<T> &getVelocitiesX() const { return this->velocitiesX; }

  /**
   * @brief Getter function for outy, the matrix containing velocities in
   * y-Direction; return a reference to outy
   */
  const RowMajMat<T> &getVelocitiesY() const { return this->velocitiesY; }

  /**
   * @brief Simulation of hydraulic charge either until convergence,
   * or for a number of selected timesteps. Calculation of Darcy-velocities.
   */
  void run() {
    // if iterations < 1 calculate hydraulic charge until steady state is
    // reached

    if constexpr (hyd_mode == STEADY_STATE) {
      // Calculate largest possible timestep, depending on K and gridsize
      const T deltaColSquare = grid.getDeltaCol() * grid.getDeltaCol();
      const T deltaRowSquare = grid.getDeltaRow() * grid.getDeltaRow();
      const T minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);

      const T maxK =
          std::max(grid.getAlphaX().maxCoeff(), grid.getAlphaY().maxCoeff());

      setTimestep(minDeltaSquare / (4 * maxK));

      RowMajMat<T> oldConcentrations;
      do {
        oldConcentrations = grid.getConcentrations();
        (void)calculate_hydraulic_flow();
      } while (!checkConvergance(oldConcentrations, grid.getConcentrations()));
    } else {
      if (timestep == -1) {
        throw_invalid_argument("Timestep is not set");
      }
      for (int i = 0; i < innerIterations; i++) {
        (void)calculate_hydraulic_flow();
      }
    }

    (void)computeFluidVelocities();
  };

private:
  /**
   * @brief Calculate the new hydraulic charge using FTCS
   */
  void calculate_hydraulic_flow() {
    if constexpr (hyd_resolve == EXPLICIT) {
      FTCS_2D(this->grid, this->bc, this->timestep, this->numThreads);
    } else {
      BTCS_2D(this->grid, this->bc, this->timestep, ThomasAlgorithm);
    }
  };

  /**
   * @brief checks if the matrix of Hydraulic Heads has converged
   *
   * @return bool true if for all corresponding cells of the matrices,
   * containing old and new Charge values, the relative error is below the
   * selected Epsilon
   */
  bool checkConvergance(const RowMajMat<T> &oldHeads,
                        const RowMajMat<T> &newHeads) {
    const auto abs_err = (oldHeads - newHeads).cwiseAbs();
    const auto rel_err = abs_err.cwiseQuotient(newHeads);

    return rel_err.maxCoeff() < epsilon;
  }

  /**
   * @brief Update the matrices containing Darcy velocities in x and y
   * directions
   */
  void computeFluidVelocities() {
    const std::size_t rows = grid.getRow();
    const std::size_t cols = grid.getCol();
    const T dx = grid.getDeltaRow();
    const T dy = grid.getDeltaCol();
    const RowMajMat<T> &hydraulicCharges = grid.getConcentrations();

    const RowMajMat<T> &permKX = grid.getAlphaX();
    const RowMajMat<T> &permKY = grid.getAlphaY();

    // calculate velocities in x-direction
    for (std::size_t i_rows = 0; i_rows < rows; i_rows++) {
      const auto bc_left = bc.getBoundaryElement(BC_SIDE_LEFT, i_rows);
      switch (bc_left.getType()) {
      case BC_TYPE_CLOSED: {
        velocitiesX(i_rows, 0) = 0;
        break;
      }
      case BC_TYPE_CONSTANT: {
        velocitiesX(i_rows, 0) =
            -permKX(i_rows, 0) *
            (hydraulicCharges(i_rows, 0) - bc_left.getValue()) / (dx / 2);
        break;
      }
      }

      const auto bc_right = bc.getBoundaryElement(BC_SIDE_RIGHT, i_rows);
      switch (bc_right.getType()) {
      case BC_TYPE_CLOSED: {
        velocitiesX(i_rows, cols) = 0;
        break;
      }
      case BC_TYPE_CONSTANT: {
        velocitiesX(i_rows, cols) =
            -permKX(i_rows, cols - 1) *
            (bc_right.getValue() - hydraulicCharges(i_rows, cols - 1)) /
            (dx / 2);
        break;
      }
      }
    }

// main loop for calculating velocities in x-direction for inner cells
#pragma omp parallel for num_threads(numThreads)
    for (int i = 0; i < rows; i++) {
      for (int j = 1; j < cols; j++) {
        velocitiesX(i, j) =
            -permKX(i, j - 1) *
            (hydraulicCharges(i, j) - hydraulicCharges(i, j - 1)) / dx;
      }
    }

    // calculate velocities in y-direction
    for (std::size_t i_cols = 0; i_cols < cols; i_cols++) {
      const auto bc_top = bc.getBoundaryElement(BC_SIDE_TOP, i_cols);
      switch (bc_top.getType()) {
      case BC_TYPE_CLOSED: {
        velocitiesY(0, i_cols) = 0;
        break;
      }
      case BC_TYPE_CONSTANT: {
        velocitiesY(0, i_cols) =
            -permKY(0, i_cols) *
            (hydraulicCharges(0, i_cols) - bc_top.getValue()) / (dy / 2);
        break;
      }
      }

      const auto bc_bottom = bc.getBoundaryElement(BC_SIDE_BOTTOM, i_cols);
      switch (bc_bottom.getType()) {
      case BC_TYPE_CLOSED: {
        velocitiesY(rows, i_cols) = 0;
        break;
      }
      case BC_TYPE_CONSTANT: {
        velocitiesY(rows, i_cols) =
            -permKY(rows - 1, i_cols) *
            (bc_bottom.getValue() - hydraulicCharges(rows - 1, i_cols)) /
            (dy / 2);
        break;
      }
      }
    }

// main loop for calculating velocities in y-direction for inner cells
#pragma omp parallel for num_threads(numThreads)
    for (int i = 1; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        velocitiesY(i, j) =
            -permKY(i, j) *
            (hydraulicCharges(i, j) - hydraulicCharges(i - 1, j)) / dy;
      }
    }
  };
};
} // namespace tug
