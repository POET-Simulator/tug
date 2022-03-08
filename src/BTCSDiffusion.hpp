#ifndef BTCSDIFFUSION_H_
#define BTCSDIFFUSION_H_

#include "BoundaryCondition.hpp"

#include <Eigen/Sparse>
#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <vector>

namespace Diffusion {
/*!
 * Class implementing a solution for a 1/2/3D diffusion equation using backward
 * euler.
 */
class BTCSDiffusion {

public:
  /*!
   * Creates a diffusion module.
   *
   * @param dim Number of dimensions. Should not be greater than 3 and not less
   * than 1.
   */
  BTCSDiffusion(unsigned int dim);

  /*!
   * Define the grid in x direction.
   *
   * @param domain_size Size of the domain in x direction.
   * @param n_grid_cells Number of grid cells in x direction the domain is
   * divided to.
   */
  void setXDimensions(double domain_size, unsigned int n_grid_cells);

  /*!
   * Define the grid in y direction.
   *
   * Throws an error if the module wasn't initialized at least as a 2D model.
   *
   * @param domain_size Size of the domain in y direction.
   * @param n_grid_cells Number of grid cells in y direction the domain is
   * divided to.
   */
  void setYDimensions(double domain_size, unsigned int n_grid_cells);

  /*!
   * Define the grid in z direction.
   *
   * Throws an error if the module wasn't initialized at least as a 3D model.
   *
   * @param domain_size Size of the domain in z direction.
   * @param n_grid_cells Number of grid cells in z direction the domain is
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
   * @param c Pointer to continious memory describing the current concentration state of each grid cell.
   * @param alpha Pointer to memory area of diffusion coefficients for each grid element.
   * @param bc Pointer to memory setting boundary conidition of each grid cell.
   */
  void simulate(double *c, double *alpha, Diffusion::boundary_condition *bc);

  /*!
   * Set the timestep of the simulation
   *
   * @param time_step Time step (in seconds ???)
   */
  void setTimestep(double time_step);

private:
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
      DMatrixRowMajor;
  typedef Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>
      DVectorRowMajor;
  typedef Eigen::Matrix<Diffusion::boundary_condition, Eigen::Dynamic,
                        Eigen::Dynamic, Eigen::RowMajor>
      BCMatrixRowMajor;
  typedef Eigen::Matrix<Diffusion::boundary_condition, 1, Eigen::Dynamic,
                        Eigen::RowMajor>
      BCVectorRowMajor;

  void simulate_base(DVectorRowMajor &c, const BCVectorRowMajor &bc,
                     const DVectorRowMajor &alpha, double dx, double time_step,
                     int size, const DVectorRowMajor &t0_c);

  void simulate1D(Eigen::Map<DVectorRowMajor> &c,
                  Eigen::Map<const DVectorRowMajor> &alpha,
                  Eigen::Map<const BCVectorRowMajor> &bc);

  void simulate2D(Eigen::Map<DMatrixRowMajor> &c,
                  Eigen::Map<const DMatrixRowMajor> &alpha,
                  Eigen::Map<const BCMatrixRowMajor> &bc);

  auto calc_t0_c(const DMatrixRowMajor &c, const DMatrixRowMajor &alpha,
                 const BCMatrixRowMajor &bc, double time_step, double dx)
      -> DMatrixRowMajor;

  inline void fillMatrixFromRow(const DVectorRowMajor &alpha,
                                const BCVectorRowMajor &bc, int size, double dx,
                                double time_step);
  inline void fillVectorFromRow(const DVectorRowMajor &c,
                                const DVectorRowMajor &alpha,
                                const BCVectorRowMajor &bc,
                                const DVectorRowMajor &t0_c, int size,
                                double dx, double time_step);
  void simulate3D(std::vector<double> &c);

  inline void reserveMemory(int size, int max_count_per_line);
  inline static auto getBCFromFlux(Diffusion::boundary_condition bc,
                                   double neighbor_c, double neighbor_alpha)
      -> double;

  void solveLES();
  void updateInternals();

  Eigen::SparseMatrix<double> A_matrix;
  Eigen::VectorXd b_vector;
  Eigen::VectorXd x_vector;

  double time_step;
  unsigned int grid_dim;

  std::vector<unsigned int> grid_cells;
  std::vector<double> domain_size;
  std::vector<double> deltas;
};
} // namespace Diffusion
#endif // BTCSDIFFUSION_H_
