/**
 * @file Velocities.hpp
 * @brief API of Velocities class, holding information for a simulation of
 * Hydraulic Charge and Darcy-Velocities.
 *
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>

#include <tug/Boundary.hpp>
#include <tug/Core/BaseSimulation.hpp>
#include <tug/Core/Matrix.hpp>
#include <tug/Core/Numeric/BTCS.hpp>
#include <tug/Core/Numeric/FTCS.hpp>
#include <tug/Core/Numeric/SimulationInput.hpp>
#include <tug/Core/TugUtils.hpp>
#include <tug/Diffusion/Diffusion.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs() 1
#endif

namespace tug {

enum class HYDRAULIC_MODE { TRANSIENT, STEADY_STATE };
enum class HYDRAULIC_RESOLVE { EXPLICIT, IMPLICIT };

template <class T, HYDRAULIC_MODE hyd_mode = HYDRAULIC_MODE::STEADY_STATE,
          HYDRAULIC_RESOLVE hyd_resolve = HYDRAULIC_RESOLVE::IMPLICIT>
class Velocities : public BaseSimulationGrid<T> {
private:
  int innerIterations{1};
  T timestep{-1};
  T epsilon{1E-5};
  int numThreads{omp_get_num_procs()};
  bool steady{false};
  bool customVelocities{false};

  RowMajMat<T> velocitiesX;
  RowMajMat<T> velocitiesY;

  RowMajMat<T> permKX;
  RowMajMat<T> permKY;

  static constexpr T DEFAULT_PERM_K = 1E-8;

  void init_velocities() {
    this->velocitiesX = RowMajMat<T>::Zero(this->rows(), this->cols() + 1);
    this->permKX =
        RowMajMat<T>::Constant(this->rows(), this->cols(), DEFAULT_PERM_K);
    if (this->getDim() == 2) {
      this->velocitiesY = RowMajMat<T>::Zero(this->rows() + 1, this->cols());
      this->permKY =
          RowMajMat<T>::Constant(this->rows(), this->cols(), DEFAULT_PERM_K);
    }
  }

public:
  /**
   * @brief Construct a new Velocities object from an Eigen RowMajMat (1D or
   * 2D).
   *
   * @param origin Reference to the hydraulic heads matrix.
   */
  Velocities(RowMajMat<T> &origin) : BaseSimulationGrid<T>(origin) {
    init_velocities();
  }

  /**
   * @brief Construct a new 2D Velocities object from raw memory pointer.
   *
   * @param data Pointer to hydraulic heads array.
   * @param rows Number of rows.
   * @param cols Number of columns.
   */
  Velocities(T *data, std::size_t rows, std::size_t cols)
      : BaseSimulationGrid<T>(data, rows, cols) {
    init_velocities();
  }

  /**
   * @brief Construct a new 1D Velocities object from raw memory pointer.
   *
   * @param data Pointer to hydraulic heads array.
   * @param length Length of the 1D grid.
   */
  Velocities(T *data, std::size_t length)
      : BaseSimulationGrid<T>(data, 1, length) {
    init_velocities();
  }

  /**
   * @brief Set the epsilon value, the relative error allowed for convergence.
   *
   * @param epsilon The new epsilon value (0 <= epsilon < 1).
   */
  void setEpsilon(T epsilon) {
    tug_assert(0 <= epsilon && epsilon < 1,
               "Relative Error epsilon must be between 0 and 1");

    this->epsilon = epsilon;
  }

  /**
   * @brief Get the permeability matrix in x-direction.
   *
   * @return RowMajMat<T>& Reference to the permKX matrix.
   */
  RowMajMat<T> &getPermKX() { return permKX; }

  const RowMajMat<T> &getPermKX() const { return permKX; }

  /**
   * @brief Get the permeability matrix in y-direction (2D only).
   *
   * @return RowMajMat<T>& Reference to the permKY matrix.
   */
  RowMajMat<T> &getPermKY() {
    tug_assert(
        this->getDim() == 2,
        "Grid is not two dimensional, there is no domain in y-direction!");

    return permKY;
  }

  const RowMajMat<T> &getPermKY() const {
    tug_assert(
        this->getDim() == 2,
        "Grid is not two dimensional, there is no domain in y-direction!");

    return permKY;
  }

  /**
   * @brief Set the permeability matrix in x-direction.
   *
   * @param alphaX The new permKX matrix.
   */
  void setPermKX(const RowMajMat<T> &alphaX) {
    tug_assert(alphaX.rows() == this->rows() && alphaX.cols() == this->cols(),
               "Permeability matrix permKX must match grid dimensions!");
    this->permKX = alphaX;
    this->steady = false;
    this->customVelocities = false;
  }

  /**
   * @brief Set the permeability matrix in y-direction (2D only).
   *
   * @param alphaY The new permKY matrix.
   */
  void setPermKY(const RowMajMat<T> &alphaY) {
    tug_assert(
        this->getDim() == 2,
        "Grid is not two dimensional, there is no domain in y-direction!");
    tug_assert(alphaY.rows() == this->rows() && alphaY.cols() == this->cols(),
               "Permeability matrix permKY must match grid dimensions!");

    this->permKY = alphaY;
    this->steady = false;
    this->customVelocities = false;
  }

  /**
   * @brief Directly set the Darcy velocities in x-direction without solving the
   * hydraulic head equation.
   *
   * @param veloX Matrix of dimensions rows() x (cols() + 1).
   */
  void setVelocitiesX(const RowMajMat<T> &veloX) {
    tug_assert(veloX.rows() == this->rows() &&
                   veloX.cols() == this->cols() + 1,
               "VelocitiesX matrix must have dimensions rows x (cols + 1)!");
    this->velocitiesX = veloX;
    this->customVelocities = true;
    this->steady = true;
  }

  /**
   * @brief Directly set the Darcy velocities in y-direction (2D only) without
   * solving the hydraulic head equation.
   *
   * @param veloY Matrix of dimensions (rows() + 1) x cols().
   */
  void setVelocitiesY(const RowMajMat<T> &veloY) {
    tug_assert(
        this->getDim() == 2,
        "Grid is not two dimensional, there is no velocity in y-direction!");
    tug_assert(veloY.rows() == this->rows() + 1 &&
                   veloY.cols() == this->cols(),
               "VelocitiesY matrix must have dimensions (rows + 1) x cols!");
    this->velocitiesY = veloY;
    this->customVelocities = true;
    this->steady = true;
  }

  bool isSteady() const { return steady; }

  /**
   * @brief Set the timestep per iteration.
   *
   * @param timestep Timestep per iteration (must be > 0).
   */
  void setTimestep(T timestep) override {
    if (timestep <= 0) {
      throw std::invalid_argument("Timestep must be greater than zero");
    }
    this->timestep = timestep;

    if constexpr (hyd_resolve == HYDRAULIC_RESOLVE::EXPLICIT) {
      const T deltaColSquare = this->deltaCol() * this->deltaCol();
      T minDeltaSquare = deltaColSquare;
      T maxK = this->permKX.maxCoeff();

      if (this->getDim() == 2) {
        const T deltaRowSquare = this->deltaRow() * this->deltaRow();
        minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);
        maxK = std::max(maxK, this->permKY.maxCoeff());
      }

      const T cfl = minDeltaSquare / (4 * maxK);
      if (timestep > cfl) {
        this->innerIterations = static_cast<int>(std::ceil(timestep / cfl));
        this->timestep = timestep / static_cast<T>(this->innerIterations);
        if (this->getOutputConsole() > CONSOLE_OUTPUT::OFF) {
          std::cerr << "Warning :: Timestep was adjusted, because of stability "
                       "conditions. Time duration was approximately preserved "
                       "by adjusting internal number of iterations."
                    << std::endl;
          std::cout << "FTCS_" << static_cast<int>(this->getDim())
                    << "D :: Required " << this->innerIterations
                    << " inner iterations with dt=" << this->timestep
                    << std::endl;
        }
      } else {
        this->innerIterations = 1;
      }
    } else {
      this->innerIterations = 1;
    }
  }

  /**
   * @brief Get the currently configured timestep.
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
      throw std::invalid_argument(
          "Number of threads exceeds the number of processor cores (" +
          std::to_string(maxThreadNumber) + ") or is less than 1.");
    }
  }

  /**
   * @brief Getter function for the matrix containing velocities in x-direction.
   *
   * @return const RowMajMat<T>& Reference to velocitiesX.
   */
  const RowMajMat<T> &getVelocitiesX() const { return this->velocitiesX; }

  /**
   * @brief Getter function for the matrix containing velocities in y-direction
   * (2D only).
   *
   * @return const RowMajMat<T>& Reference to velocitiesY.
   */
  const RowMajMat<T> &getVelocitiesY() const {
    tug_assert(
        this->getDim() == 2,
        "Grid is not two dimensional, there is no velocity in y-direction!");
    return this->velocitiesY;
  }

  /**
   * @brief Simulation of hydraulic charge either until convergence,
   * or for a number of selected timesteps. Calculation of Darcy-velocities.
   */
  void run() override {
    if (this->customVelocities) {
      return;
    }

    this->applyInnerBoundaries();

    SimulationInput<T> input = {
        .concentrations = this->getConcentrationMatrix(),
        .alphaX = this->permKX,
        .alphaY = this->permKY,
        .boundaries = this->getBoundaryConditions(),
        .dim = this->getDim(),
        .timestep = this->timestep,
        .rowMax = this->rows(),
        .colMax = this->cols(),
        .deltaRow =
            (this->getDim() == 2) ? this->deltaRow() : static_cast<T>(0),
        .deltaCol = this->deltaCol()};

    if constexpr (hyd_mode == HYDRAULIC_MODE::STEADY_STATE) {
      if (steady) {
        return;
      }

      const T deltaColSquare = this->deltaCol() * this->deltaCol();
      T minDeltaSquare = deltaColSquare;
      T maxK = this->permKX.maxCoeff();

      if (this->getDim() == 2) {
        const T deltaRowSquare = this->deltaRow() * this->deltaRow();
        minDeltaSquare = std::min(deltaColSquare, deltaRowSquare);
        maxK = std::max(maxK, this->permKY.maxCoeff());
      }

      const T cflTimestep = minDeltaSquare / (4 * maxK);
      if (this->getOutputConsole() > CONSOLE_OUTPUT::OFF) {
        std::cout << "CFL (hydHead) timestep: " << cflTimestep << std::endl;
      }

      setTimestep(cflTimestep);
      input.timestep = this->timestep;

      RowMajMat<T> oldConcentrations;
      do {
        oldConcentrations = this->getConcentrationMatrix();
        calculate_hydraulic_flow(input);
      } while (!checkConvergence(oldConcentrations));

      steady = true;

    } else {
      if (timestep <= 0) {
        throw_invalid_argument("Timestep is not set");
      }
      for (int i = 0; i < innerIterations; i++) {
        calculate_hydraulic_flow(input);
      }
    }

    computeFluidVelocities();
  }

private:
  /**
   * @brief Calculate the new hydraulic charge using FTCS or BTCS.
   */
  void calculate_hydraulic_flow(SimulationInput<T> &sim_in) {
    if constexpr (hyd_resolve == HYDRAULIC_RESOLVE::EXPLICIT) {
      FTCS(sim_in, numThreads);
    } else {
      BTCS_Thomas(sim_in, numThreads);
    }
  }

  /**
   * @brief Checks if the matrix of Hydraulic Heads has converged.
   *
   * @return bool True if for all corresponding cells of the matrices,
   * containing old and new Charge values, the relative error is below the
   * selected Epsilon.
   */
  bool checkConvergence(const RowMajMat<T> &oldHeads) const {
    const auto &newHeads = this->getConcentrationMatrix();
    const auto abs_err = (oldHeads - newHeads).cwiseAbs();
    const auto denom = newHeads.cwiseAbs().cwiseMax(epsilon);
    const auto rel_err = abs_err.cwiseQuotient(denom);

    return rel_err.maxCoeff() < epsilon;
  }

  /**
   * @brief Update the matrices containing Darcy velocities in x and y
   * directions.
   */
  void computeFluidVelocities() {
    const std::size_t rows = this->rows();
    const std::size_t cols = this->cols();
    const T dx = this->deltaCol();
    const RowMajMat<T> &hydraulicCharges = this->getConcentrationMatrix();

    const RowMajMat<T> &permKX = this->permKX;
    const Boundary<T> &bc = this->getBoundaryConditions();

    // calculate velocities in x-direction at domain boundaries
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

    // main loop for calculating velocities in x-direction for inner cell edges
#pragma omp parallel for num_threads(numThreads)
    for (std::size_t i = 0; i < rows; i++) {
      for (std::size_t j = 1; j < cols; j++) {
        const T k_inter = calcAlphaIntercell(permKX(i, j - 1), permKX(i, j));
        velocitiesX(i, j) =
            -k_inter *
            (hydraulicCharges(i, j) - hydraulicCharges(i, j - 1)) / dx;
      }
    }

    if (this->getDim() == 2) {
      const T dy = this->deltaRow();
      const RowMajMat<T> &permKY = this->permKY;

      // calculate velocities in y-direction at domain boundaries
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

      // main loop for calculating velocities in y-direction for inner cell
      // edges
#pragma omp parallel for num_threads(numThreads)
      for (std::size_t i = 1; i < rows; i++) {
        for (std::size_t j = 0; j < cols; j++) {
          const T k_inter = calcAlphaIntercell(permKY(i - 1, j), permKY(i, j));
          velocitiesY(i, j) =
              -k_inter *
              (hydraulicCharges(i, j) - hydraulicCharges(i - 1, j)) / dy;
        }
      }
    }
  }
};
} // namespace tug
