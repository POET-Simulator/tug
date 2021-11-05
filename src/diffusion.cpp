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

      tripletList.push_back(T(A_line, i * x + j, (1 + 2 * sx + 2 * sy)));
      tripletList.push_back(T(A_line, (i - 1) * x + j, sy));
      tripletList.push_back(T(A_line, (i + 1) * x + j, sy));
      tripletList.push_back(T(A_line, i * x + (j + 1), sx));
      tripletList.push_back(T(A_line, i * x + (j - 1), sx));

      b[A_line] = -c[i*x+j];
      A_line++;
    }
  }

  std::cout << b << std::endl;

  Eigen::SparseMatrix<double> A(size, (x * y) - 4);
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
