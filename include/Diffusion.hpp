#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "BoundaryCondition.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <bits/stdint-uintn.h>
#include <vector>

namespace tug {
namespace diffusion {

/**
 * Defines grid dimensions and boundary conditions.
 */
typedef struct {
  uint32_t grid_cells[3];
  double domain_size[3];
  bc::BoundaryCondition *bc;
} TugGrid;

/**
 * Besides containing the grid structure it holds also information about the
 * desired time step to simulate and which solver to use.
 */
typedef struct {
  double time_step;
  Eigen::VectorXd (*solver)(Eigen::SparseMatrix<double>, Eigen::VectorXd);
  TugGrid grid;
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
