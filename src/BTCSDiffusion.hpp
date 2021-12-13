#ifndef BTCSDIFFUSION_H_
#define BTCSDIFFUSION_H_

#include <Eigen/Sparse>
#include <Eigen/src/SparseCore/SparseMatrixBase.h>
#include <tuple>
#include <vector>

/*!
 * Datatype to fill the sparse matrix which is used to solve the equation
 * system.
 */
typedef Eigen::Triplet<double> T;

/*!
 * Defines both types of boundary condition as a datatype.
 */
typedef int bctype;

/*!
 * A boundary condition consists of two features. A type and the according
 * value. Here we can differentiate between:
 *
 * - Neumann boundary conditon: type BC_NEUMANN with the value defining the
 * gradient
 * - Dirichlet boundary condition: type BC_DIRICHLET with the actual value of
 * the boundary condition
 */
typedef std::vector<std::tuple<bctype, double>> boundary_condition;

/*!
 * Class implementing a solution for a 1/2/3D diffusion equation using backward
 * euler.
 */
class BTCSDiffusion {

public:
  /*!
   * Defines a Neumann boundary condition.
   */
  static const int BC_NEUMANN;
  /*!
   * Defines a Dirichlet boundary condition.
   */
  static const int BC_DIRICHLET;

  /*!
   * Create 1D-diffusion module.
   *
   * @param x Count of cells in x direction.
   */
  BTCSDiffusion(int x);

  /*!
   * Currently not implemented: Create 2D-diffusion module.
   *
   * @param x Count of cells in x direction.
   * @param y Count of cells in y direction.
   */
  BTCSDiffusion(int x, int y);

  /*!
   * Currently not implemented: Create 3D-diffusion module.
   *
   * @param x Count of cells in x direction.
   * @param y Count of cells in y direction.
   * @param z Count of cells in z direction.
   */
  BTCSDiffusion(int x, int y, int z);

  /*!
   * With given ghost zones simulate diffusion. Only 1D allowed at this moment.
   *
   * @param c Vector describing the concentration of one solution of the grid as
   * continious memory (Row-wise).
   * @param alpha Vector of diffusioncoefficients for each grid element.
   */
  void simulate(std::vector<double> &c, std::vector<double> &alpha);

  /*!
   * Set the timestep of the simulation
   *
   * @param time_step Time step (in seconds ???)
   */
  void setTimestep(double time_step);

  /*!
   * Set the boundary condition of the given grid. This is done by defining an
   * index (exact order still to be determined), the type of the boundary
   * condition and the according value.
   *
   * @param index Index of the boundary condition vector.
   * @param val Value of the boundary condition (gradient for Neumann, exact
   * value for Dirichlet).
   * @param Type of the grid cell.
   */
  void setBoundaryCondition(int index, double val, bctype type);

private:
  void simulate1D(std::vector<double> &c, double bc_left, double bc_right,
                  std::vector<double> &alpha);
  void simulate2D(std::vector<double> &c);
  void simulate3D(std::vector<double> &c);

  double getBCFromTuple(int index, double nearest_value);

  boundary_condition bc;

  Eigen::SparseMatrix<double> A_matrix;
  Eigen::VectorXd b_vector;
  Eigen::VectorXd x_vector;

  double time_step;

  int grid_dim;
  int dim_x;
  int dim_y;
  int dim_z;
};

#endif // BTCSDIFFUSION_H_
