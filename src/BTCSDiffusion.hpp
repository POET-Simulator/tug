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
  unsigned int getXGridCellsN();
  /*!
   * Returns the number of grid cells in y direction.
   */
  unsigned int getYGridCellsN();
  /*!
   * Returns the number of grid cells in z direction.
   */
  unsigned int getZGridCellsN();

  /*!
   * Returns the domain size in x direction.
   */
  unsigned int getXDomainSize();
  /*!
   * Returns the domain size in y direction.
   */
  unsigned int getYDomainSize();
  /*!
   * Returns the domain size in z direction.
   */
  unsigned int getZDomainSize();

  /*!
   * With given ghost zones simulate diffusion. Only 1D allowed at this moment.
   *
   * @param c Vector describing the concentration of one solution of the grid as
   * continious memory (row major).
   * @param alpha Vector of diffusion coefficients for each grid element.
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
   * @param index Index of the grid cell the boundary condition is applied to.
   * @param type Type of the boundary condition. Must be constant, closed or
   * flux.
   * @param value For constant boundary conditions this value is set
   * during solving. For flux value refers to a gradient of change for this grid
   * cell. For closed this value has no effect since a gradient of 0 is used.
   */
  void setBoundaryCondition(int index, bctype type, double value);

private:
  typedef struct boundary_condition {
    bctype type;
    double value;
  } boundary_condition;
  typedef Eigen::Triplet<double> T;

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
