#include "diffusion/BTCSBoundaryCondition.hpp"
#include <diffusion/BTCSDiffusion.hpp>
#include <iomanip>
#include <iostream> // for std
#include <vector>   // for vector
using namespace std;

using namespace Diffusion;

int main(int argc, char *argv[]) {

  // dimension of grid
  int dim = 2;

  int n = 5;
  int m = 5;

  // create input + diffusion coefficients for each grid cell
  std::vector<double> alpha(n * m, 1e-6);
  std::vector<double> field(n * m, 0);
  BTCSBoundaryCondition bc(n, m);

  // create instance of diffusion module
  BTCSDiffusion diffu(dim);

  diffu.setXDimensions(n, n);
  diffu.setYDimensions(m, m);

  // set inlet to higher concentration without setting bc
  field[12] = 1;

  // set timestep for simulation to 1 second
  diffu.setTimestep(1);

  cout << setprecision(12);

  for (int t = 0; t < 1000; t++) {
    diffu.simulate(field.data(), alpha.data(), bc);

    cout << "Iteration: " << t << "\n\n";

    double sum = 0;

    // iterate through rows
    for (int i = 0; i < m; i++) {
      // iterate through columns
      for (int j = 0; j < n; j++) {
        cout << left << std::setw(20) << field[i * n + j];
        sum += field[i * n + j];
      }
      cout << "\n";
    }

    cout << "sum: " << sum << "\n" << endl;
  }

  return 0;
}
