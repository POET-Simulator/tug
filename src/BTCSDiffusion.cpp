#include "BTCSDiffusion.hpp"

#include <Eigen/SparseLU>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

const int BTCSDiffusion::BC_NEUMANN = 0;
const int BTCSDiffusion::BC_DIRICHLET = 1;

BTCSDiffusion::BTCSDiffusion(int x) : dim_x(x) {
  this->grid_dim = 1;

  // per default use Neumann condition with gradient of 0 at the end of the grid
  this->bc.resize(2, std::tuple<bctype, double>(0, 0.));
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

  /*
   * Begin to solve the equation system using LU solver of Eigen.
   *
   * But first fill the A matrix and b vector.
   *
   * At this point there is some debugging output in the code.
   * TODO: remove output
   */

  A_matrix.resize(size, size);
  A_matrix.reserve(Eigen::VectorXi::Constant(size, 3));

  A_matrix.insert(0, 0) = 1;
  A_matrix.insert(size - 1, size - 1) = 1;

  b_vector[0] = bc_left;
  b_vector[size - 1] = bc_right;

  for (int i = 1; i < this->dim_x + 1; i++) {
    double sx = (alpha[i - 1] * time_step) / (dx * dx);

    A_matrix.insert(i, i) = -1. - 2. * sx;
    A_matrix.insert(i, i - 1) = sx;
    A_matrix.insert(i, i + 1) = sx;

    b_vector[i] = -c[i - 1];

  }

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
    double bc_left = getBCFromTuple(0, c[0]);
    double bc_right = getBCFromTuple(1, c[c.size() - 1]);
    // double bc_left = 5. * std::pow(10,-6);
    // double bc_right = c[this->dim_x -1];

    simulate1D(c, bc_left, bc_right, alpha);
  }
}

double BTCSDiffusion::getBCFromTuple(int index, double nearest_value) {
  double val = -1;
  int type = std::get<0>(bc[index]);

  if (type == BTCSDiffusion::BC_NEUMANN) {
    // TODO implement gradient here
    val = nearest_value;
  } else if (type == BTCSDiffusion::BC_DIRICHLET) {
    val = std::get<1>(bc[index]);
  } else {
    // TODO: implement error handling here. Type was set to wrong value.
  }

  return val;
}

void BTCSDiffusion::setBoundaryCondition(int index, double val, bctype type) {
  std::get<0>(bc[index]) = type;
  std::get<1>(bc[index]) = val;
}
