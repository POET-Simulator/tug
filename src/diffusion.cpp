#include "diffusion.hpp"

#include <Eigen/SparseLU>
#include <Eigen/src/Core/Matrix.h>

void BTCS2D(int x, int y, std::vector<double> &c, std::vector<double> &alpha,
            double timestep) {

  double dx = 1. / x;
  double dy = 1. / y;

  int size = x * y;

  Eigen::VectorXd d = Eigen::VectorXd::Constant(size, 0);
  std::vector<T> tripletList;
  tripletList.reserve(((x - 1) * (y - 1) * 5));

  for (int i = 1; i < y - 1; i++) {
    for (int j = 1; j < x - 1; j++) {
      double sx = (alpha[i * y + j] * timestep) / (dx * dx);
      double sy = (alpha[i * y + j] * timestep) / (dy * dy);

      tripletList.push_back(T((i * x) + j, (i * x) + j, (1 + 2 * sx + 2 * sy)));
      tripletList.push_back(T((i * x) + j, ((i * x) + j) + x, sy));
      tripletList.push_back(T((i * x) + j, ((i * x) + j) - x, sy));
      tripletList.push_back(T((i * x) + j, ((i * x) + j) + 1, sx));
      tripletList.push_back(T((i * x) + j, ((i * x) + j) - 1, sx));
    }
  }
}
