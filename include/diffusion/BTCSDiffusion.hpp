#ifndef BTCSDIFFUSION_H_
#define BTCSDIFFUSION_H_

#include "grid/BoundaryCondition.hpp"

#include <Eigen/Sparse>
#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <vector>

namespace tug {
namespace diffusion {

/*!
 * Class implementing a solution for a 1/2/3D diffusion equation using backward
 * euler.
 */
class BTCSDiffusion {

public:
  /*!
   * Creates a diffusion module.
   *
   * \param dim Number of dimensions. Should not be greater than 3 and not less
   * than 1.
   */
  BTCSDiffusion(unsigned int dim);

  /*!
   * Define the grid in x direction.
   *
   * \param domain_size Size of the domain in x direction.
   * \param n_grid_cells Number of grid cells in x direction the domain is
   * divided to.
   */
  void setXDimensions(double domain_size, unsigned int n_grid_cells);

  /*!
   * Define the grid in y direction.
   *
   * Throws an error if the module wasn't initialized at least as a 2D model.
   *
   * \param domain_size Size of the domain in y direction.
   * \param n_grid_cells Number of grid cells in y direction the domain is
   * divided to.
   */
  void setYDimensions(double domain_size, unsigned int n_grid_cells);

  /*!
   * Define the grid in z direction.
   *
   * Throws an error if the module wasn't initialized at least as a 3D model.
   *
   * \param domain_size Size of the domain in z direction.
   * \param n_grid_cells Number of grid cells in z direction the domain is
   * divided to.
   */
  void setZDimensions(double domain_size, unsigned int n_grid_cells);

  /*!
   * Returns the number of grid cells in x direction.
   */
  auto getXGridCellsN() -> unsigned int;
  /*!
   * Returns the number of grid cells in y direction.
   */
  auto getYGridCellsN() -> unsigned int;
  /*!
   * Returns the number of grid cells in z direction.
   */
  auto getZGridCellsN() -> unsigned int;

  /*!
   * Returns the domain size in x direction.
   */
  auto getXDomainSize() -> double;
  /*!
   * Returns the domain size in y direction.
   */
  auto getYDomainSize() -> double;
  /*!
   * Returns the domain size in z direction.
   */
  auto getZDomainSize() -> double;

  /*!
   * With given ghost zones simulate diffusion. Only 1D allowed at this moment.
   *
   * \param c Pointer to continious memory describing the current concentration
   * state of each grid cell.
   * \param alpha Pointer to memory area of diffusion coefficients for each grid
   * element.
   * \param bc Instance of boundary condition class with.
   *
   * \return Time in seconds [s] used to simulate one iteration.
   */
  auto simulate(double *c, double *alpha,
                const tug::boundary_condition::BoundaryCondition &bc)
      -> double;

  /*!
   * Set the timestep of the simulation
   *
   * \param time_step Time step (in seconds ???)
   */
  void setTimestep(double time_step);

private:
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      DMatrixRowMajor;
  typedef Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>
      DVectorRowMajor;

  static void simulate_base(DVectorRowMajor &c,
                            const tug::boundary_condition::bc_tuple &bc_ghosts,
                            const tug::boundary_condition::bc_vec &bc_inner,
                            const DVectorRowMajor &alpha, double dx,
                            double time_step, int size,
                            const DVectorRowMajor &d_ortho);

  void simulate1D(Eigen::Map<DVectorRowMajor> &c,
                  Eigen::Map<const DVectorRowMajor> &alpha,
                  const tug::boundary_condition::BoundaryCondition &bc);

  void simulate2D(Eigen::Map<DMatrixRowMajor> &c,
                  Eigen::Map<const DMatrixRowMajor> &alpha,
                  const tug::boundary_condition::BoundaryCondition &bc);

  static auto
  calc_d_ortho(const DMatrixRowMajor &c, const DMatrixRowMajor &alpha,
               const tug::boundary_condition::BoundaryCondition &bc,
               bool transposed, double time_step, double dx) -> DMatrixRowMajor;

  static void fillMatrixFromRow(Eigen::SparseMatrix<double> &A_matrix,
                                const DVectorRowMajor &alpha,
                                const tug::boundary_condition::bc_vec &bc_inner,
                                int size, double dx, double time_step);

  static void
  fillVectorFromRow(Eigen::VectorXd &b_vector, const DVectorRowMajor &c,
                    const DVectorRowMajor &alpha,
                    const tug::boundary_condition::bc_tuple &bc_ghosts,
                    const tug::boundary_condition::bc_vec &bc_inner,
                    const DVectorRowMajor &d_ortho, int size, double dx,
                    double time_step);
  void simulate3D(std::vector<double> &c);

  inline static auto
  getBCFromFlux(tug::boundary_condition::boundary_condition bc,
                double neighbor_c, double neighbor_alpha) -> double;

  void updateInternals();

  double time_step;
  unsigned int grid_dim;

  std::vector<unsigned int> grid_cells;
  std::vector<double> domain_size;
  std::vector<double> deltas;
};
} // namespace diffusion
} // namespace tug
#endif // BTCSDIFFUSION_H_
