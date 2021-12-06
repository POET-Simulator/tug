#include "BTCSDiffusion.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

  int x = 20;

  std::vector<double> alpha(x, 1 * pow(10, -1));
  std::vector<double> input(x, 1 * std::pow(10, -6));
  std::vector<double> bc_left, bc_right;

  bc_left.push_back(5. * std::pow(10, -6));
  bc_right.push_back(-1);

  BTCSDiffusion diffu(x);

  diffu.setBoundaryCondition(0, 5. * std::pow(10, -6), BTCSDiffusion::BC_DIRICHLET);
  diffu.setTimestep(1.);

  // diffu.setBoundaryCondition(bc_left, BTCSDiffusion::LEFT);
  // we don't need this since Neumann condition with gradient of 0 is set per
  // default
  // diffu.setBoundaryCondition(bc_right, BTCSDiffusion::RIGHT);

  for (int i = 0; i < 100; i++) {
    diffu.simulate(input, alpha);
  }

  return 0;
}
