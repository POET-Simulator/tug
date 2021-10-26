#include "diffusion.hpp"

#include <Eigen/SparseLU>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/OrderingMethods/Ordering.h>
#include <Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include <Eigen/src/SparseCore/SparseMap.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseMatrixBase.h>
#include <Eigen/src/SparseLU/SparseLU.h>
#include<Eigen/SparseCholesky>
#include<Eigen/SparseQR>
#include <Eigen/src/SparseQR/SparseQR.h>
#include <iostream>

void BTCS2D(int x, int y, std::vector<double> &c, std::vector<double> &alpha,
            double timestep) {

  double dx = 1. / x;
  double dy = 1. / y;

  int size = (x * y) - 4;

  int local_x = x - 2;

  Eigen::VectorXd b = Eigen::VectorXd::Constant(size, 0);
  Eigen::VectorXd x_out(size);
  std::vector<T> tripletList;
  tripletList.reserve(size * 5);

  for (int i = x - 1 ; i < 2*x - 3 ; i++) {
    double sx = (alpha[i + 2] * timestep) / (dx * dx);
    double sy = (alpha[i + 2] * timestep) / (dy * dy);


    tripletList.push_back(T(i, i, (1 + 2 * sx + 2 * sy)));
    tripletList.push_back(T(i, i - (x - 1), sy));
    tripletList.push_back(T(i, i + x, sy));
    tripletList.push_back(T(i, i + 1, sx));
    tripletList.push_back(T(i, i - 1, sx));

    b[i] = -c[i+2];
  }

  for (int i = 2*x - 1; i < (y-2)*x-3; i++) {
    double sx = (alpha[i + 2] * timestep) / (dx * dx);
    double sy = (alpha[i + 2] * timestep) / (dy * dy);

    tripletList.push_back(T(i, i, (1 + 2 * sx + 2 * sy)));
    tripletList.push_back(T(i, i - x, sy));
    tripletList.push_back(T(i, i + x, sy));
    tripletList.push_back(T(i, i + 1, sx));
    tripletList.push_back(T(i, i - 1, sx));

    b[i] = -c[i+2];
  }

  for (int i = (y-2)*x-1; i < (y-1)*x-3; i++) {
    double sx = (alpha[i + 2] * timestep) / (dx * dx);
    double sy = (alpha[i + 2] * timestep) / (dy * dy);

    tripletList.push_back(T(i, i, (1 + 2 * sx + 2 * sy)));
    tripletList.push_back(T(i, i - x, sy));
    tripletList.push_back(T(i, i + (x-1), sy));
    tripletList.push_back(T(i, i + 1, sx));
    tripletList.push_back(T(i, i - 1, sx));

    b[i] = -c[i+2];
  }

  // for (int i = 0; i < (size-local_x); i++) {
  //   int current = local_x + i;
  //   double sx = (alpha[current] * timestep) / (dx * dx);
  //   double sy = (alpha[current] * timestep) / (dy * dy);

  //   tripletList.push_back(T(i, current, (1 + 2 * sx + 2 * sy)));
  //   tripletList.push_back(T(i, current + local_x, sy));
  //   tripletList.push_back(T(i, current - local_x, sy));
  //   tripletList.push_back(T(i, current + 1, sx));
  //   tripletList.push_back(T(i, current - 1, sx));

  //   std::cout << current << std::endl;

  //   b[i] = -c[x+i];
  // }

  // for (int i = 0; i < y; i++) {
  //   for (int j = 0; j < x; j++) {
  //     double sx = (alpha[i * y + j] * timestep) / (dx * dx);
  //     double sy = (alpha[i * y + j] * timestep) / (dy * dy);

  //     tripletList.push_back(T((i * x) + j, (i * x) + j, (1 + 2 * sx + 2 *
  //     sy))); tripletList.push_back(T((i * x) + j, ((i * x) + j) + x, sy));
  //     tripletList.push_back(T((i * x) + j, ((i * x) + j) - x, sy));
  //     tripletList.push_back(T((i * x) + j, ((i * x) + j) + 1, sx));
  //     tripletList.push_back(T((i * x) + j, ((i * x) + j) - 1, sx));

  //     b[(i * x) + j] = -c[(i * x) + j];
  //   }
  // }

  std::cout << b << std::endl;

  Eigen::SparseMatrix<double> A(size, (x * y) - 4);
  A.setFromTriplets(tripletList.begin(), tripletList.end());

Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

  // Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
  //     solver;
  solver.analyzePattern(A);

  solver.factorize(A);

  std::cout << A << std::endl;
  std::cout << solver.lastErrorMessage() << std::endl;

  x_out = solver.solve(b);

  std::cout << x_out << std::endl;
}
