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
  std::vector<double> alpha(n * m, 1e-1);
  std::vector<double> field(n * m, 1e-6);
  BTCSBoundaryCondition bc(n, m);

  // create instance of diffusion module
  BTCSDiffusion diffu(dim);

  diffu.setXDimensions(1, n);
  diffu.setYDimensions(1, m);

  boundary_condition input = {Diffusion::BC_TYPE_CONSTANT, 5e-6};

  bc.setSide(BC_SIDE_LEFT, input);

  // for (int i = 1; i <= n; i++) {
  //   bc[(n + 2) * i] = {Diffusion::BC_CONSTANT, 5e-6};
  // }
  // set timestep for simulation to 1 second
  diffu.setTimestep(1.);

  cout << setprecision(12);

  for (int t = 0; t < 10; t++) {
    diffu.simulate(field.data(), alpha.data(), bc);

    cout << "Iteration: " << t << "\n\n";

    // iterate through rows
    for (int i = 0; i < m; i++) {
      // iterate through columns
      for (int j = 0; j < n; j++) {
        cout << left << std::setw(20) << field[i * n + j];
      }
      cout << "\n";
    }

    cout << "\n" << endl;
  }

  return 0;
}
