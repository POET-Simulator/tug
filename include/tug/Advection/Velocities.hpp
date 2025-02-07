/**
 * @file Velocities.hpp
 * @brief API of Velocities class, holding information for a simulation of
 * Hydraulic Charge and Darcy-Velocities. Holds a predifined Grid object and
 * Boundary object.
 *
 */

#pragma once

#include "tug/Core/Numeric/SimulationInput.hpp"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <string>

#include <tug/Boundary.hpp>
#include <tug/Core/BaseSimulation.hpp>
#include <tug/Core/Matrix.hpp>
#include <tug/Core/Numeric/BTCS.hpp>
#include <tug/Core/Numeric/FTCS.hpp>
#include <tug/Core/TugUtils.hpp>
#include <tug/Diffusion/Diffusion.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs() 1
#endif

using namespace Eigen;
namespace tug {

enum class HYDRAULIC_MODE { TRANSIENT, STEADY_STATE };
enum class HYDRAULIC_RESOLVE { EXPLICIT, IMPLICIT };

template <class T, HYDRAULIC_MODE hyd_mode, HYDRAULIC_RESOLVE hyd_resolve>
class Velocities : public BaseSimulationGrid<T> {
private:
  int innerIterations{1};
  T timestep{-1};
  T epsilon{1E-5};
  int numThreads{omp_get_num_procs()};
  bool steady{false};

  RowMajMat<T> velocitiesX;
  RowMajMat<T> velocitiesY;

  RowMajMat<T> permKX;
  RowMajMat<T> permKY;

public:
  Velocities(RowMajMat<T> &origin)
      : BaseSimulationGrid<T>(origin),
        velocitiesX(origin.rows(), origin.cols() + 1),
        velocitiesY(origin.rows() + 1, origin.cols()),
        permKX(origin.rows(), origin.cols()),
        permKY(origin.rows(), origin.cols()) {};

  Velocities(T *data, std::size_t rows, std::size_t cols)
      : BaseSimulationGrid<T>(data, rows, cols), velocitiesX(rows, cols + 1),
        velocitiesY(rows + 1, cols), permKX(rows, cols), permKY(rows, cols) {};

  // Velocities(T *data, std::size_t length)
  //     : BaseSimulationGrid<T>(data, 1, length), velocitiesX(1, length + 1),
  //       alphaX(1, length) {};

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
   * @brief Get the alphaX matrix.
   *
   * @return RowMajMat<T>& Reference to the alphaX matrix.
   */
  RowMajMat<T> &getPermKX() { return permKX; }

  /**
   * @brief Get the alphaY matrix.
   *
   * @return RowMajMat<T>& Reference to the alphaY matrix.
   */
  RowMajMat<T> &getPermKY() {
    tug_assert(
        this->getDim(),
        "Grid is not two dimensional, there is no domain in y-direction!");

    return permKY;
  }

  /**
   * @brief Set the alphaX matrix.
   *
   * @param alphaX The new alphaX matrix.
   */
  void setPermKX(const RowMajMat<T> &alphaX) { this->permKX = alphaX; }

  /**
   * @brief Set the alphaY matrix.
   *
   * @param alphaY The new alphaY matrix.
   */
  void setPermKY(const RowMajMat<T> &alphaY) {
    tug_assert(
        this->getDim(),
        "Grid is not two dimensional, there is no domain in y-direction!");

    this->permKY = alphaY;
  }

  bool isSteady() const { return steady; }

  /**
   * @brief Set the timestep per iteration
   *
   * @param timestep timestep per iteration
   */
  void setTimestep(T timestep) override {
    if (timestep <= 0) {
      throw std::invalid_argument("Timestep must be greater than zero");
    }
    this->timestep = timestep;
    const T deltaColSquare = this->deltaCol() * this->deltaCol();
    const T deltaRowSquare = this->deltaRow() * this->deltaRow();
    const T minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);

    const T maxK = std::max(this->permKX.maxCoeff(), this->permKY.maxCoeff());

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
  void run() override {
    // if iterations < 1 calculate hydraulic charge until steady state is
    // reached

    SimulationInput<T> input = {.concentrations =
                                    this->getConcentrationMatrix(),
                                .alphaX = this->getPermKX(),
                                .alphaY = this->getPermKY(),
                                .boundaries = this->getBoundaryConditions(),
                                .dim = this->getDim(),
                                .timestep = this->timestep,
                                .rowMax = this->rows(),
                                .colMax = this->cols(),
                                .deltaRow = this->deltaRow(),
                                .deltaCol = this->deltaCol()};

    if constexpr (hyd_mode == HYDRAULIC_MODE::STEADY_STATE) {
      if (steady) {
        return;
      }

      const T deltaColSquare = this->deltaCol() * this->deltaCol();
      const T deltaRowSquare = this->deltaRow() * this->deltaRow();
      const T minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);

      const T maxK = std::max(this->permKX.maxCoeff(), this->permKY.maxCoeff());
      // Calculate largest possible timestep, depending on K and gridsize

      setTimestep(minDeltaSquare / (4 * maxK));

      input.timestep = this->timestep;

      RowMajMat<T> oldConcentrations;
      do {
        oldConcentrations = this->getConcentrationMatrix();
        (void)calculate_hydraulic_flow(input);
      } while (!checkConvergance(oldConcentrations));

      steady = true;

    } else {
      if (timestep == -1) {
        throw_invalid_argument("Timestep is not set");
      }
      for (int i = 0; i < innerIterations; i++) {
        (void)calculate_hydraulic_flow(input);
      }
    }

    (void)computeFluidVelocities();
  };

private:
  /**
   * @brief Calculate the new hydraulic charge using FTCS
   */
  void calculate_hydraulic_flow(SimulationInput<T> &sim_in) {
    if constexpr (hyd_resolve == HYDRAULIC_RESOLVE::EXPLICIT) {
      FTCS_2D(sim_in, numThreads);
    } else {
      BTCS_2D(sim_in, ThomasAlgorithm, numThreads);
    }
  };

  /**
   * @brief checks if the matrix of Hydraulic Heads has converged
   *
   * @return bool true if for all corresponding cells of the matrices,
   * containing old and new Charge values, the relative error is below the
   * selected Epsilon
   */
  bool checkConvergance(const RowMajMat<T> &oldHeads) {
    const auto abs_err = (oldHeads - this->getConcentrationMatrix()).cwiseAbs();
    const auto rel_err = abs_err.cwiseQuotient(this->getConcentrationMatrix());

    return rel_err.maxCoeff() < epsilon;
  }

  /**
   * @brief Update the matrices containing Darcy velocities in x and y
   * directions
   */
  void computeFluidVelocities() {
    const std::size_t rows = this->rows();
    const std::size_t cols = this->cols();
    const T dx = this->deltaCol();
    const T dy = this->deltaRow();
    const RowMajMat<T> &hydraulicCharges = this->getConcentrationMatrix();

    const RowMajMat<T> &permKX = this->permKX;
    const RowMajMat<T> &permKY = this->permKY;

    const Boundary<T> &bc = this->getBoundaryConditions();

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
            -this->permKX(i_rows, 0) *
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
