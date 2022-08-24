#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "BoundaryCondition.hpp"
#include "Solver.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <bits/stdint-uintn.h>
#include <vector>

namespace tug {
namespace diffusion {

constexpr uint8_t MAX_ARR_SIZE = 3;

/**
 * Defines grid dimensions and boundary conditions.
 */
typedef struct tug_grid_s {
  std::array<uint32_t, MAX_ARR_SIZE>
      grid_cells; /**< Count of grid cells in each of the 3 directions.*/
  std::array<double, MAX_ARR_SIZE>
      domain_size; /**< Domain sizes in each of the 3 directions.*/
  bc::BoundaryCondition *bc = NULL; /**< Boundary conditions for the grid.*/
} TugGrid;

/**
 * Besides containing the grid structure it holds also information about the
 * desired time step to simulate and which solver to use.
 */
typedef struct tug_input_s {
  double time_step; /**< Time step which should be simulated by diffusion.*/
  Eigen::VectorXd (*solver)(const Eigen::SparseMatrix<double> &,
                            const Eigen::VectorXd &) =
      tug::solver::ThomasAlgorithm; /**< Solver function to use.*/
  TugGrid grid;                     /**< Grid specification.*/

  /**
   * Set the desired time step for diffusion simulation.
   *
   * \param dt Time step in seconds.
   */
  void setTimestep(double dt) { time_step = dt; }

  /**
   * Set the count of grid cells in each dimension.
   *
   * \param x Count of grid cells in x direction.
   * \param y Count of grid cells in y direction.
   * \param z Count of grid cells in z direction.
   */
  void setGridCellN(uint32_t x, uint32_t y = 0, uint32_t z = 0) {
    grid.grid_cells[0] = x;
    grid.grid_cells[1] = y;
    grid.grid_cells[2] = z;
  }

  /**
   * Set the domain size of the grid in each direction.

   * \param Domain size in x direction.
   * \param Domain size in y direction.
   * \param Domain size in z direction.
   */
  void setDomainSize(double x, double y = 0, double z = 0) {
    grid.domain_size[0] = x;
    grid.domain_size[1] = y;
    grid.domain_size[2] = z;
  }

  /**
   * Set boundary conditions for grid instance.
   *
   * \param bc Boundary conditions to be set.
   */
  void setBoundaryCondition(bc::BoundaryCondition &bc) { grid.bc = &bc; }

  /**
   * Retrieve the set boundary condition from grid instance.
   *
   * \return Boundary condition object if boundary conditions were set,
   * otherwise NULL.
   */
  auto getBoundaryCondition() -> bc::BoundaryCondition { return *(grid.bc); }

  /**
   * Set the solver function.
   *
   * \param f_in Pointer to function which takes a sparse matrix and a vector as
   * input and returns another vector.
   */
  void
  setSolverFunction(Eigen::VectorXd (*f_in)(const Eigen::SparseMatrix<double> &,
                                            const Eigen::VectorXd &)) {
    solver = f_in;
  }
} TugInput;

/**
 * Solving 1D diffusion problems with Backward Time Centred Space scheme.
 *
 * \param input_param Object with information for the simulation e.g. grid
 * dimensions, time step ...
 *
 * \param field Pointer to continious memory holding the concentrations for each
 * grid cell of the grid. It doesn't matter if stored in row (most likely) or
 * column major.
 *
 * \param alpha Pointer to continious memory holding the alpha for each grid
 * cell of the grid. (NOTE: only constant alphas are supported currently)
 *
 * \return Runtime in seconds
 */
auto BTCS_1D(const TugInput &input_param, double *field, const double *alpha)
    -> double;

/**
 * Solving 2D diffusion problems with Alternating-direction implicit method.
 *
 * \param input_param Object with information for the simulation e.g. grid
 * dimensions, time step ...
 *
 * \param field Pointer to continious memory holding the concentrations for each
 * grid cell of the grid. It doesn't matter if stored in row (most likely) or
 * column major.
 *
 * \param alpha Pointer to continious memory holding the alpha for each grid
 * cell of the grid. (NOTE: only constant alphas are supported currently)
 *
 * \return Runtime in seconds
 */
auto ADI_2D(const TugInput &input_param, double *field, const double *alpha)
    -> double;

} // namespace diffusion
} // namespace tug

#endif // DIFFUSION_H_
