#include "BTCSDiffusion.hpp"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <tuple>
#include <vector>

const int BTCSDiffusion::BC_CONSTANT = 0;
const int BTCSDiffusion::BC_CLOSED = 1;
const int BTCSDiffusion::BC_FLUX = 2;

BTCSDiffusion::BTCSDiffusion(unsigned int dim) : grid_dim(dim) {
  assert(dim <= 3);

  grid_cells.resize(dim, 1);
  domain_size.resize(dim, 1);
  deltas.resize(dim, 1);
}

std::vector<int> BTCSDiffusion::getNumberOfGridCells() {
  return this->grid_cells;
}
std::vector<int> BTCSDiffusion::getSpatialDiscretization() {
  return this->domain_size;
}
void BTCSDiffusion::setNumberOfGridCells(std::vector<int> &n_grid) {
  grid_cells = n_grid;
  assert(grid_cells.size() == grid_dim);
  updateInternals();
}
void BTCSDiffusion::setSpatialDiscretization(std::vector<int> &s_grid) {
  domain_size = s_grid;
  assert(domain_size.size() == grid_dim);
  updateInternals();
}

void BTCSDiffusion::updateInternals() {
  for (int i = 0; i < grid_dim; i++) {
    deltas[i] = (double)domain_size[i] / grid_cells[i];
  }

  switch (grid_dim) {
  case 1:
    bc.resize(2, {BTCSDiffusion::BC_CLOSED, 0});
    break;
  case 2:
    bc.resize(2 * grid_cells[0] + 2 * grid_cells[1],
              {BTCSDiffusion::BC_CLOSED, 0});
    break;
  case 3:
    // TODO
    break;
  }
}
// BTCSDiffusion::BTCSDiffusion(int x) : n_x(x) {
//   this->grid_dim = 1;
//   this->dx = 1. / (x - 1);

//   // per default use Neumann condition with gradient of 0 at the end of the
//   grid this->bc.resize(2, std::tuple<bctype,
//   double>(BTCSDiffusion::BC_CONSTANT, 0.));
// }
// BTCSDiffusion::BTCSDiffusion(int x, int y) : n_x(x), n_y(y) {

//   // this->grid_dim = 2;

//   // this->bc.reserve(x * 2 + y * 2);
//   // // per default use Neumann condition with gradient of 0 at the end of
//   the
//   // grid std::fill(this->bc.begin(), this->bc.end(), -1);
// }
// BTCSDiffusion::BTCSDiffusion(int x, int y, int z) : n_x(x), n_y(y), n_z(z) {

//   // this->grid_dim = 3;
//   // TODO: reserve memory for boundary conditions
// }

void BTCSDiffusion::simulate1D(std::vector<double> &c, boundary_condition left,
                               boundary_condition right,
                               const std::vector<double> &alpha, double dx,
                               int size) {

  bool left_is_constant = (left.type == BTCSDiffusion::BC_CONSTANT);
  bool right_is_constant = (right.type == BTCSDiffusion::BC_CONSTANT);
  int loop_end = size + !right_is_constant;

  // we need 2 more grid cells for ghost cells
  // size = size + 2;

  int bc_offset = !left_is_constant + !right_is_constant;
  ;

  // set sizes of private and yet allocated vectors
  b_vector.resize(size + bc_offset);
  x_vector.resize(size + bc_offset);

  /*
   * Begin to solve the equation system using LU solver of Eigen.
   *
   * But first fill the A matrix and b vector.
   *
   * At this point there is some debugging output in the code.
   * TODO: remove output
   */

  A_matrix.resize(size + bc_offset, size + bc_offset);
  A_matrix.reserve(Eigen::VectorXi::Constant(size + bc_offset, 3));

  A_matrix.insert(0, 0) = 1;
  b_vector[0] = (left_is_constant ? left.value : getBCFromFlux(left, c[0], alpha[0]));

  A_matrix.insert((size + bc_offset) - 1, (size + bc_offset) - 1) = 1;
  b_vector[size + bc_offset - 1] =
    (right_is_constant ? right.value : getBCFromFlux(right, c[size - 1], alpha[size - 1]));

  // A_matrix.insert(0, 0) = 1;
  // A_matrix.insert(size + 1, size + 1) = 1;

  for (int i = 1; i < size - right_is_constant; i++) {
    double sx = (alpha[i + !(left_is_constant)] * time_step) / (dx * dx);

    A_matrix.insert(i, i) = -1. - 2. * sx;
    A_matrix.insert(i, i - 1) = sx;
    A_matrix.insert(i, i + 1) = sx;

    b_vector[i] = -c[i + !(left_is_constant)];
  }

std::cout << b_vector << "\n" << A_matrix << std::endl;

  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A_matrix);

  solver.factorize(A_matrix);

  std::cout << solver.lastErrorMessage() << std::endl;

  // x_vector = solver.solve(b_vector);

  std::cout << std::setprecision(10) << x_vector << std::endl << std::endl;

  for (int i = 0; i < c.size(); i++) {
    c[i] = x_vector[i + !left_is_constant];
  }
}

void BTCSDiffusion::setTimestep(double time_step) {
  this->time_step = time_step;
}

void BTCSDiffusion::simulate(std::vector<double> &c,
                             const std::vector<double> &alpha) {
  if (this->grid_dim == 1) {
    // double bc_left = getBCFromTuple(0, c[0], alpha[0]);
    // double bc_right =
    //     getBCFromTuple(1, c[c.size() - 1], alpha[alpha.size() - 1]);

    simulate1D(c, bc[0], bc[1], alpha, this->deltas[0], this->grid_cells[0]);
  }
}

inline double BTCSDiffusion::getBCFromFlux(boundary_condition bc, double neighbor_c,
                                    double neighbor_alpha) {

  double val;

  if (bc.type == BTCSDiffusion::BC_CLOSED) {
    val = neighbor_c;
  } else if (bc.type == BTCSDiffusion::BC_FLUX) {
    //TODO
    // val = bc[index].value;
  } else {
    // TODO: implement error handling here. Type was set to wrong value.
  }

  return val;
}

void BTCSDiffusion::setBoundaryCondition(int index, double val, bctype type) {

  bc[index].type = type;
  bc[index].value = val;

  // std::get<0>(bc[index]) = type;
  // std::get<1>(bc[index]) = val;
}
