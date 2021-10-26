#include "diffusion.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

  int x = 10;
  int y = 10;

  std::vector<double> alpha(x * y, 1.5 * pow(10, -9));
  std::vector<double> input(x * y, 0);
  input[x + 1] = 5 * std::pow(10, 6);

  BTCS2D(x, y, input, alpha, 10.);

  return 0;
}
