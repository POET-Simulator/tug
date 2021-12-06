#include "BTCSDiffusion.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/OrderingMethods/Ordering.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <tuple>

const int BTCSDiffusion::BC_NEUMANN = 0;
const int BTCSDiffusion::BC_DIRICHLET = 1;

BTCSDiffusion::BTCSDiffusion(int x) : dim_x(x) {
  this->grid_dim = 1;

  // per default use Neumann condition with gradient of 0 at the end of the grid
  this->bc.resize(2, std::tuple<int, double>(0, 0.));
}
BTCSDiffusion::BTCSDiffusion(int x, int y) : dim_x(x), dim_y(y) {

  // this->grid_dim = 2;

  // this->bc.reserve(x * 2 + y * 2);
  // // per default use Neumann condition with gradient of 0 at the end of the
  // grid std::fill(this->bc.begin(), this->bc.end(), -1);
}
BTCSDiffusion::BTCSDiffusion(int x, int y, int z)
    : dim_x(x), dim_y(y), dim_z(z) {

  // this->grid_dim = 3;
  // TODO: reserve memory for boundary conditions
}

void BTCSDiffusion::simulate1D(std::vector<double> &c, double bc_left,
                               double bc_right, std::vector<double> &alpha) {
  // calculate dx
  double dx = 1. / (this->dim_x - 1);

  // calculate size needed for A matrix and b,x vectors
  int size = this->dim_x + 2;

  // set sizes of private and yet allocated vectors
  b_vector.resize(size);
  x_vector.resize(size);

  // Eigen::VectorXd b = Eigen::VectorXd::Constant(size, 0);
  // Eigen::VectorXd x_out(size);

  /*
   * Initalization of matrix A
   * This is done by triplets. See:
   * https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
   */

  // std::vector<T> tripletList;
  // tripletList.reserve(c.size() * 3 + bc.size());

  // int A_line = 0;

  // // For all concentrations create one row in matrix A
  // for (int i = 1; i < this->dim_x + 1; i++) {
  //   double sx = (alpha[i - 1] * timestep) / (dx * dx);

  //   tripletList.push_back(T(A_line, i, (-1. - 2. * sx)));

  //   tripletList.push_back(T(A_line, i - 1, sx));
  //   tripletList.push_back(T(A_line, i + 1, sx));

  //   b[A_line] = -c[i - 1];
  //   A_line++;
  // }

  // // append left and right boundary conditions/ghost zones
  // tripletList.push_back(T(A_line, 0, 1));

  // // if value is -1 apply Neumann condition with given gradient
  // // TODO: set specific gradient
  // if (bc[0] == -1)
  //   b[A_line] = c[0];
  // // else apply given Dirichlet condition
  // else
  //   b[A_line] = this->bc[0];

  // A_line++;
  // tripletList.push_back(T(A_line, size - 1, 1));
  // // b[A_line] = bc[1];
  // if (bc[1] == -1)
  //   b[A_line] = c[c.size() - 1];
  // else
  //   b[A_line] = this->bc[1];

  /*
   * Begin to solve the equation system
   *
   * At this point there is some debugging output in the code.
   * TODO: remove output
   */

  A_matrix.resize(size, size);
  A_matrix.reserve(Eigen::VectorXi::Constant(size, 3));

  A_matrix.insert(0, 0) = bc_left;
  A_matrix.insert(size - 1, size - 1) = bc_right;

  for (int i = 1; i < this->dim_x + 1; i++) {
    double sx = (alpha[i - 1] * time_step) / (dx * dx);

    A_matrix.insert(i, i) = -1. - 2. * sx;
    A_matrix.insert(i, i - 1) = sx;
    A_matrix.insert(i, i + 1) = sx;

    b_vector[i] = -c[i - 1];

    // tripletList.push_back(T(A_line, i, (-1. - 2. * sx)));

    // tripletList.push_back(T(A_line, i - 1, sx));
    // tripletList.push_back(T(A_line, i + 1, sx));

    // b[A_line] = -c[i - 1];
    // A_line++;
  }

  // Eigen::SparseMatrix<double> A(size, size);
  // A.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A_matrix);

  solver.factorize(A_matrix);

  std::cout << solver.lastErrorMessage() << std::endl;

  x_vector = solver.solve(b_vector);

  std::cout << std::setprecision(10) << x_vector << std::endl << std::endl;

  for (int i = 0; i < c.size(); i++) {
    c[i] = x_vector[i + 1];
  }
}

void BTCSDiffusion::setTimestep(double time_step) {
  this->time_step = time_step;
}

void BTCSDiffusion::simulate(std::vector<double> &c,
                             std::vector<double> &alpha) {
  if (this->grid_dim == 1) {
    double bc_left = getBCFromTuple(0);
    double bc_right = getBCFromTuple(1);

    simulate1D(c, bc_left, bc_right, alpha);
  }
}

double BTCSDiffusion::getBCFromTuple(int index) {
  double val = std::get<1>(bc[index]);

  return val;
}

void BTCSDiffusion::setBoundaryCondition(int index, double val, int type) {
  std::get<0>(bc[index]) = val;
  std::get<1>(bc[index]) = type;
}
