#include "diffusion.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

  int x = 4;
  int y = 4;

  std::vector<double> alpha(x * y, 1 * pow(10, -9));
  std::vector<double> input(x * y, 1 * std::pow(10,-6));
  input[x + 1] = 5 * std::pow(10, -6);
  // input[x + 2] = 5.5556554 * std::pow(10, -6);
  // input[x + 3] = 5.234564213 * std::pow(10, -6);

  BTCS2D(x, y, input, alpha, 1.);

  return 0;
}
