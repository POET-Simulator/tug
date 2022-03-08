#include "BTCSDiffusion.hpp" // for BTCSDiffusion, BTCSDiffusion::BC_DIRICHLET
#include "BoundaryCondition.hpp"
#include <algorithm>         // for copy, max
#include <cmath>
#include <iomanip>
#include <iostream> // for std
#include <vector>   // for vector
using namespace std;

using namespace Diffusion;

int main(int argc, char *argv[]) {

  // dimension of grid
  int dim = 2;

  int n = 501;
  int m = 501;

  // create input + diffusion coefficients for each grid cell
  std::vector<double> alpha(n * m, 1 * pow(10, -1));
  std::vector<double> field(n * m, 0.);
  std::vector<boundary_condition> bc(n*m, {0,0});

  field[125500] = 1;

  // create instance of diffusion module
  BTCSDiffusion diffu(dim);

  diffu.setXDimensions(1., n);
  diffu.setYDimensions(1., m);

  // set the boundary condition for the left ghost cell to dirichlet
  //diffu.setBoundaryCondition(250, 250, BTCSDiffusion::BC_CONSTANT, 1);
  // for (int d=0; d<5;d++){
  //     diffu.setBoundaryCondition(d, 0, BC_CONSTANT, .1);
  // }
  // diffu.setBoundaryCondition(1, 1, BTCSDiffusion::BC_CONSTANT, .1);
  // diffu.setBoundaryCondition(1, 1, BTCSDiffusion::BC_CONSTANT, .1);
  
  // set timestep for simulation to 1 second
  diffu.setTimestep(1.);

  cout << setprecision(7);

  // First we output the initial state
  cout << 0;
    
  for (int i=0; i < m*n; i++) {
      cout << "," << field[i];
  }
  cout << endl;

  // Now we simulate and output 8 steps à 1 sec
  for (int t = 1; t < 6; t++) {
    diffu.simulate(field.data(), alpha.data(), bc.data());
    cout << t;
    
    for (int i=0; i < m*n; i++) {
        cout << "," << field[i];
    }
    cout << endl;
  }

  return 0;
}
