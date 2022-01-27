#ifndef BTCSDIFFUSION_H_
#define BTCSDIFFUSION_H_

#include <Eigen/Sparse>
#include <tuple>
#include <type_traits>
#include <vector>

/*!
 * Defines both types of boundary condition as a datatype.
 */
typedef int bctype;

/*!
 * Class implementing a solution for a 1/2/3D diffusion equation using backward
 * euler.
 */
class BTCSDiffusion {

public:
  /*!
   * Defines a constant/Dirichlet boundary condition.
   */
  static const int BC_CONSTANT;

  /*!
   * Defines a closed/Neumann boundary condition.
   */
  static const int BC_CLOSED;

  /*!
   * Defines a flux/Cauchy boundary condition.
   */
  static const int BC_FLUX;

  /*!
   * A boundary condition consists of two features. A type and the according
   * value. Here we can differentiate between:
   *
   * - Neumann boundary conditon: type BC_NEUMANN with the value defining the
   * gradient
   * - Dirichlet boundary condition: type BC_DIRICHLET with the actual value of
   * the boundary condition
   */
  typedef struct boundary_condition {
    bctype type;
    double value;
  } boundary_condition;

  /*!
   * A boundary condition consists of two features. A type and the according
   * value. Here we can differentiate between:
   *
   * - Neumann boundary conditon: type BC_NEUMANN with the value defining the
   * gradient
   * - Dirichlet boundary condition: type BC_DIRICHLET with the actual value of
   * the boundary condition
   */
  // typedef std::vector<std::tuple<bctype, double>> boundary_condition;

  /*!
   * Datatype to fill the sparse matrix which is used to solve the equation
   * system.
   */
  typedef Eigen::Triplet<double> T;

  /*!
   * Create 1D-diffusion module.
   *
   * @param x Count of cells in x direction.
   */
  BTCSDiffusion(unsigned int dim);

  void setXDimensions(unsigned int domain_size, unsigned int n_grid_cells);
  void setYDimensions(unsigned int domain_size, unsigned int n_grid_cells);
  void setZDimensions(unsigned int domain_size, unsigned int n_grid_cells);

  unsigned int getXGridCellsN();
  unsigned int getYGridCellsN();
  unsigned int getZGridCellsN();
  unsigned int getXDomainSize();
  unsigned int getYDomainSize();
  unsigned int getZDomainSize();
  // /*!
  //  * Currently not implemented: Create 2D-diffusion module.
  //  *
  //  * @param x Count of cells in x direction.
  //  * @param y Count of cells in y direction.
  //  */
  // explicit BTCSDiffusion(int x, int y);

  // /*!
  //  * Currently not implemented: Create 3D-diffusion module.
  //  *
  //  * @param x Count of cells in x direction.
  //  * @param y Count of cells in y direction.
  //  * @param z Count of cells in z direction.
  //  */
  // explicit BTCSDiffusion(int x, int y, int z);

  /*!
   * With given ghost zones simulate diffusion. Only 1D allowed at this moment.
   *
   * @param c Vector describing the concentration of one solution of the grid as
   * continious memory (Row-wise).
   * @param alpha Vector of diffusioncoefficients for each grid element.
   */
  void simulate(std::vector<double> &c, const std::vector<double> &alpha);

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
  void simulate1D(std::vector<double> &c, boundary_condition left,
                  boundary_condition right, const std::vector<double> &alpha,
                  double dx, int size);
  void simulate2D(std::vector<double> &c);
  void simulate3D(std::vector<double> &c);

  inline double getBCFromFlux(boundary_condition bc, double nearest_value,
                              double neighbor_alpha);

  void updateInternals();

  std::vector<boundary_condition> bc;

  Eigen::SparseMatrix<double> A_matrix;
  Eigen::VectorXd b_vector;
  Eigen::VectorXd x_vector;

  double time_step;

  int grid_dim;
  std::vector<unsigned int> grid_cells;
  std::vector<unsigned int> domain_size;
  std::vector<double> deltas;
};

#endif // BTCSDIFFUSION_H_
