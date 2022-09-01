#include "grid/BoundaryCondition.hpp"
#include <Diffusion.hpp>
#include <iostream>
#include <vector>

#define SIZE 5

using namespace tug;

int main() {
  boundary_condition::BoundaryCondition bc(SIZE);

  diffusion::TugDiffuInput in = {1., {SIZE, 0, 0}, {1, 0, 0}, 0};

  std::vector<double> field1(SIZE, 0);
  std::vector<double> field2(SIZE, 0);
  std::vector<double> alpha(SIZE, 1e-5);

  bc(boundary_condition::BC_SIDE_LEFT) = {boundary_condition::BC_TYPE_CONSTANT,
                                          1};

  for (int j = 0; j < 10; j++) {
    double time1 =
        diffusion::BTCS_Thomas_1D(field1.data(), alpha.data(), bc, in);
    double time2 =
        diffusion::BTCS_EigenLU_1D(field2.data(), alpha.data(), bc, in);

    for (int i = 0; i < SIZE; i++) {
      std::cout << field1[i] << "\t";
      std::cout << field2[i] << "\n";
    }

    std::cout << "Time Thomas = " << time1 << "s \tTime EigenLU = " << time2
              << " s" << std::endl;
  }
}
