#include "BTCSDiffusion.hpp"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/OrderingMethods/Ordering.h>
#include <Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include <Eigen/src/SparseCore/SparseMap.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseMatrixBase.h>
#include <Eigen/src/SparseLU/SparseLU.h>
#include <Eigen/src/SparseQR/SparseQR.h>

#include <iomanip>
#include <iostream>

const BCSide BTCSDiffusion::LEFT = 0;
const BCSide BTCSDiffusion::RIGHT = 1;

BTCSDiffusion::BTCSDiffusion(int x) : dim_x(x) {
  this->grid_dim = 1;
  this->bc.reserve(2);
}
BTCSDiffusion::BTCSDiffusion(int x, int y) : dim_x(x), dim_y(y) {

  this->grid_dim = 2;
  this->bc.reserve(x * 2 + y * 2);
}
BTCSDiffusion::BTCSDiffusion(int x, int y, int z)
    : dim_x(x), dim_y(y), dim_z(z) {

  this->grid_dim = 3;
  //TODO: reserve memory for boundary conditions
}

void BTCSDiffusion::setBoundaryCondition(std::vector<double> input,
                                         BCSide side) {
  if (this->grid_dim == 1) {
    bc[side] = input[0];
  }
}
void BTCSDiffusion::simulate(std::vector<double> &c, std::vector<double> &alpha,
                             double timestep) {
  double dx = 1. / this->dim_x;

  int size = this->dim_x + 2;

  Eigen::VectorXd b = Eigen::VectorXd::Constant(size, 0);
  Eigen::VectorXd x_out(size);
  std::vector<T> tripletList;
  tripletList.reserve(c.size() * 3 + bc.size());

  int A_line = 0;

  for (int i = 1; i < this->dim_x + 1; i++) {
    double sx = (alpha[i - 1] * timestep) / (dx * dx);

    tripletList.push_back(T(A_line, i, (-1. - 2. * sx)));

    tripletList.push_back(T(A_line, i - 1, sx));
    tripletList.push_back(T(A_line, i + 1, sx));

    b[A_line] = -c[i - 1];
    A_line++;
  }

  tripletList.push_back(T(A_line, 0, 1));
  if (bc[0] == -1)
    b[A_line] = c[0];
  else
    b[A_line] = this->bc[0];

  A_line++;
  tripletList.push_back(T(A_line, size - 1, 1));
  // b[A_line] = bc[1];
  if (bc[1] == -1)
    b[A_line] = c[c.size() - 1];
  else
    b[A_line] = this->bc[1];

  Eigen::SparseMatrix<double> A(size, size);
  A.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;

  // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
  //     solver;
  solver.analyzePattern(A);

  solver.factorize(A);

  std::cout << solver.lastErrorMessage() << std::endl;

  x_out = solver.solve(b);

  std::cout << std::setprecision(10) << x_out << std::endl << std::endl;

  for (int i = 0; i < c.size(); i++) {
    c[i] = x_out[i + 1];
  }
}
