#include "diffusion.hpp"

#include <Eigen/SparseLU>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <Eigen/src/OrderingMethods/Ordering.h>
#include <Eigen/src/SparseCore/SparseMap.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseMatrixBase.h>
#include <Eigen/src/SparseLU/SparseLU.h>
#include <iostream>

void BTCS2D(int x, int y, std::vector<double> &c, std::vector<double> &alpha,
            double timestep) {

  double dx = 1. / x;
  double dy = 1. / y;

  int size = (x - 2) * (y - 2);

  int local_x = x - 2;

  Eigen::VectorXd b = Eigen::VectorXd::Constant(size, 0);
  Eigen::VectorXd x_out(size);
  std::vector<T> tripletList;
  tripletList.reserve(size * 5);

  for (int i = 0; i < (size-local_x); i++) {
    int current = local_x + i;
    double sx = (alpha[current] * timestep) / (dx * dx);
    double sy = (alpha[current] * timestep) / (dy * dy);

    tripletList.push_back(T(i, current, (1 + 2 * sx + 2 * sy)));
    tripletList.push_back(T(i, current + local_x, sy));
    tripletList.push_back(T(i, current - local_x, sy));
    tripletList.push_back(T(i, current + 1, sx));
    tripletList.push_back(T(i, current - 1, sx));

    std::cout << current << std::endl;

    b[i] = -c[x+i];
  }

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

  Eigen::SparseMatrix<double> A(size, (x*y)-4);
  A.setFromTriplets(tripletList.begin(), tripletList.end());

  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A);

  solver.factorize(A);

  std::cout << A << std::endl;

  x_out = solver.solve(b);
}
