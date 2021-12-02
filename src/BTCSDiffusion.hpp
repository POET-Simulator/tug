#ifndef BTCSDIFFUSION_H_
#define BTCSDIFFUSION_H_

#include <Eigen/Sparse>
#include <vector>

/*!
 * Type defining the side of given boundary condition.
 */
typedef int BCSide;

/*!
 * Datatype to fill the sparse matrix which is used to solve the equation
 * system.
 */
typedef Eigen::Triplet<double> T;

/*!
 * Class implementing a solution for a 1/2/3D diffusion equation using backward
 * euler.
 */
class BTCSDiffusion {

public:
  /*!
   * Set left boundary condition.
   */
  static const BCSide LEFT;

  /*!
   * Set right boundary condition.
   */
  static const BCSide RIGHT;

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
   * Sets internal boundary condition at the end of the grid/ghost zones.
   * Currently only implemented for 1D diffusion.
   *
   * @param input Vector containing all the values to initialize the ghost
   * zones.
   * @param side Sets the side of the boundary condition. See BCSide for more
   * information.
   */
  void setBoundaryCondition(std::vector<double> input, BCSide side);

  /*!
   * With given ghost zones simulate diffusion. Only 1D allowed at this moment.
   *
   * @param c Vector describing the concentration of one solution of the grid as
   * continious memory (Row-wise).
   * @param alpha Vector of diffusioncoefficients for each grid element.
   * @param timestep Time (in seconds ?) to simulate.
   */
  void simulate(std::vector<double> &c, std::vector<double> &alpha,
                double timestep);

private:
  std::vector<double> bc;
  int grid_dim;
  int dim_x;
  int dim_y;
  int dim_z;
};

#endif // BTCSDIFFUSION_H_
