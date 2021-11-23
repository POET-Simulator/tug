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
  // input[x + 2] = 5.5556554 * std::pow(10, -6);
  // input[x + 3] = 5.234564213 * std::pow(10, -6);

  BTCSDiffusion diffu(x);

  diffu.setBoundaryCondition(bc_left, BTCSDiffusion::LEFT);
  diffu.setBoundaryCondition(bc_right, BTCSDiffusion::RIGHT);

  for (int i = 0; i < 100; i++) {
    diffu.simulate(input, alpha, 1.);
    // BTCS1D(x, input, alpha, 1., bc);
  }

  return 0;
}
