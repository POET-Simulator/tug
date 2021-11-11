#include "diffusion.hpp"

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
#include <iostream>
#include <iomanip>
#include <ostream>

void BTCS1D(int x, std::vector<double> &c, std::vector<double> &alpha,
            double timestep, std::vector<double> &bc) {
  double dx = 1. / x;

  int size = x + 2;

  Eigen::VectorXd b = Eigen::VectorXd::Constant(size, 0);
  Eigen::VectorXd x_out(size);
  std::vector<T> tripletList;
  tripletList.reserve(c.size() * 3 + bc.size());

  int A_line = 0;

  for (int i = 1; i < x + 1; i++) {
    double sx = (alpha[i-1] * timestep) / (dx * dx);

    tripletList.push_back(T(A_line, i, (-1. - 2. * sx)));

    tripletList.push_back(T(A_line, i - 1, sx));
    tripletList.push_back(T(A_line, i + 1, sx));

    b[A_line] = -c[i-1];
    A_line++;
  }

  tripletList.push_back(T(A_line, 0, 1));
  b[A_line] = bc[0];

  A_line++;
  tripletList.push_back(T(A_line, size-1, 1));
  b[A_line] = bc[1];

  // std::cout << b << std::endl;

  Eigen::SparseMatrix<double> A(size, size);
  A.setFromTriplets(tripletList.begin(), tripletList.end());

  // std::cout << A << std::endl;
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;

  // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
  //     solver;
  solver.analyzePattern(A);

  solver.factorize(A);

  std::cout << solver.lastErrorMessage() << std::endl;

  x_out = solver.solve(b);

  std::cout << std::setprecision(10) << x_out << std::endl << std::endl;

  for (int i=0; i < c.size(); i++) {
    c[i] = x_out[i+1];
  }
}

void BTCS2D(int x, int y, std::vector<double> &c, std::vector<double> &alpha,
            double timestep) {

  double dx = 1. / x;
  double dy = 1. / y;

  int size = (x - 2) * (y - 2);

  Eigen::VectorXd b = Eigen::VectorXd::Constant(size, 0);
  Eigen::VectorXd x_out(x * y);
  std::vector<T> tripletList;
  tripletList.reserve(size * 5);

  int A_line = 0;

  for (int i = 1; i < y - 1; i++) {
    for (int j = 1; j < x - 1; j++) {
      double sx = (alpha[i * x + j] * timestep) / (dx * dx);
      double sy = (alpha[i * x + j] * timestep) / (dy * dy);

      tripletList.push_back(T(A_line, i * x + j, (1. + 2. * sx + 2. * sy)));

      std::cout << sx << std::endl;

      tripletList.push_back(T(A_line, (i - 1) * x + j, sy));
      tripletList.push_back(T(A_line, (i + 1) * x + j, sy));
      tripletList.push_back(T(A_line, i * x + (j + 1), sx));
      tripletList.push_back(T(A_line, i * x + (j - 1), sx));

      b[A_line] = -c[i * x + j];
      A_line++;
    }
  }

  std::cout << b << std::endl;

  Eigen::SparseMatrix<double> A(size, x * y);
  A.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;

  // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
  //     solver;
  solver.analyzePattern(A);

  solver.factorize(A);

  std::cout << A << std::endl;
  std::cout << solver.lastErrorMessage() << std::endl;

  x_out = solver.solve(b);

  std::cout << x_out << std::endl;
}
