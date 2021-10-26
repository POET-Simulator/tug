#include "diffusion.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

int main(int argc, char *argv[]) {

  int x = 5;
  int y = 5;

  std::vector<double> alpha(x * y, -1.5 * pow(10, -2));
  std::vector<double> input(x * y, 0);
  input[x + 1] = 5.55 * std::pow(10, -6);
  input[x + 2] = 5.5556554 * std::pow(10, -6);
  input[x + 3] = 5.234564213 * std::pow(10, -6);

  BTCS2D(x, y, input, alpha, 10.);

  return 0;
}
