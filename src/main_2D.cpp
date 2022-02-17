#include "BTCSDiffusion.hpp" // for BTCSDiffusion, BTCSDiffusion::BC_DIRICHLET
#include <algorithm>         // for copy, max
#include <iomanip>
#include <iostream> // for std
#include <vector>   // for vector
using namespace std;

int main(int argc, char *argv[]) {

  // dimension of grid
  int dim = 2;

  int n = 5;
  int m = 5;

  // create input + diffusion coefficients for each grid cell
  std::vector<double> alpha(n * m, 1 * pow(10, -1));
  std::vector<double> field(n * m, 1 * std::pow(10, -6));

  // create instance of diffusion module
  BTCSDiffusion diffu(dim);

  diffu.setXDimensions(1, n);
  diffu.setYDimensions(1, m);

  // set the boundary condition for the left ghost cell to dirichlet
  diffu.setBoundaryCondition(2, 2, BTCSDiffusion::BC_CONSTANT,
                             5. * std::pow(10, -6));

  // set timestep for simulation to 1 second
  diffu.setTimestep(1.);

  cout << setprecision(12);

  for (int t = 0; t < 10; t++) {
    diffu.simulate(field, alpha);

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
