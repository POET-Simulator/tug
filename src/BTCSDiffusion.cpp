#include "BTCSDiffusion.hpp"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cassert>
#include <iomanip>
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

void BTCSDiffusion::setXDimensions(unsigned int domain_size,
                                   unsigned int n_grid_cells) {
  assert(this->grid_dim > 0);
  this->domain_size[0] = domain_size;
  this->grid_cells[0] = n_grid_cells;

  updateInternals();
}

void BTCSDiffusion::setYDimensions(unsigned int domain_size,
                                   unsigned int n_grid_cells) {
  assert(this->grid_dim > 1);
  this->domain_size[1] = domain_size;
  this->grid_cells[1] = n_grid_cells;

  updateInternals();
}

void BTCSDiffusion::setZDimensions(unsigned int domain_size,
                                   unsigned int n_grid_cells) {
  assert(this->grid_dim > 2);
  this->domain_size[2] = domain_size;
  this->grid_cells[2] = n_grid_cells;

  updateInternals();
}

unsigned int BTCSDiffusion::getXGridCellsN() { return this->grid_cells[0]; }
unsigned int BTCSDiffusion::getYGridCellsN() { return this->grid_cells[1]; }
unsigned int BTCSDiffusion::getZGridCellsN() { return this->grid_cells[2]; }
unsigned int BTCSDiffusion::getXDomainSize() { return this->domain_size[0]; }
unsigned int BTCSDiffusion::getYDomainSize() { return this->domain_size[1]; }
unsigned int BTCSDiffusion::getZDomainSize() { return this->domain_size[2]; }

void BTCSDiffusion::updateInternals() {
  for (int i = 0; i < grid_dim; i++) {
    deltas[i] = (double)domain_size[i] / grid_cells[i];
  }

  int cells = 1;

  for (int i = 0; i < grid_dim; i++) {
    cells *= (grid_cells[i] + 2);
  }

  bc.resize(cells, {BTCSDiffusion::BC_CLOSED, 0});
}

void BTCSDiffusion::simulate1D(std::vector<double> &c, boundary_condition left,
                               boundary_condition right,
                               const std::vector<double> &alpha, double dx,
                               int size) {

  bool left_is_constant = (left.type == BTCSDiffusion::BC_CONSTANT);
  bool right_is_constant = (right.type == BTCSDiffusion::BC_CONSTANT);

  //The sizes for matrix and vectors of the equation system is defined by the
  //actual size of the input vector and if the system is (partially) closed.
  //Then we will need ghost nodes. So this variable will give the count of ghost
  //nodes.
  int bc_offset = !left_is_constant + !right_is_constant;
  ;

  // set sizes of private and yet allocated vectors
  b_vector.resize(size + bc_offset);
  x_vector.resize(size + bc_offset);

  /*
   * Begin to solve the equation system using LU solver of Eigen.
   *
   * But first fill the A matrix and b vector.
   */

  // Set boundary condition for ghost nodes (for closed or flux system) or outer
  // inlet nodes (constant boundary condition)
  A_matrix.resize(size + bc_offset, size + bc_offset);
  A_matrix.reserve(Eigen::VectorXi::Constant(size + bc_offset, 3));

  A_matrix.insert(0, 0) = 1;
  b_vector[0] =
      (left_is_constant ? left.value : getBCFromFlux(left, c[0], alpha[0]));

  A_matrix.insert((size + bc_offset) - 1, (size + bc_offset) - 1) = 1;
  b_vector[size + bc_offset - 1] =
      (right_is_constant ? right.value
                         : getBCFromFlux(right, c[size - 1], alpha[size - 1]));

  // Start filling the A matrix
  // =i= is used for equation system matrix and vector indexing
  // and =j= for indexing of c,alpha and bc
  for (int i = 1, j = i + !(left_is_constant); i < size - right_is_constant;
       i++, j++) {

    // if current grid cell is considered as constant boundary conditon
    if (bc[j].type == BTCSDiffusion::BC_CONSTANT) {
      A_matrix.insert(i, i) = 1;
      b_vector[i] = bc[j].value;
      continue;
    }

    double sx = (alpha[j] * time_step) / (dx * dx);

    A_matrix.insert(i, i) = -1. - 2. * sx;
    A_matrix.insert(i, i - 1) = sx;
    A_matrix.insert(i, i + 1) = sx;

    b_vector[i] = -c[j];
  }

  // start to solve
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A_matrix);

  solver.factorize(A_matrix);

  x_vector = solver.solve(b_vector);

  //fill solution back in place into =c= vector
  for (int i = 0, j = i + !left_is_constant; i < c.size(); i++, j++) {
    c[i] = x_vector[i + !left_is_constant];
  }
}

void BTCSDiffusion::setTimestep(double time_step) {
  this->time_step = time_step;
}

void BTCSDiffusion::simulate(std::vector<double> &c,
                             const std::vector<double> &alpha) {
  if (this->grid_dim == 1) {
    simulate1D(c, bc[0], bc[grid_cells[0] + 1], alpha, this->deltas[0],
               this->grid_cells[0]);
  }
}

inline double BTCSDiffusion::getBCFromFlux(boundary_condition bc,
                                           double neighbor_c,
                                           double neighbor_alpha) {

  double val;

  if (bc.type == BTCSDiffusion::BC_CLOSED) {
    val = neighbor_c;
  } else if (bc.type == BTCSDiffusion::BC_FLUX) {
    // TODO
    //  val = bc[index].value;
  } else {
    // TODO: implement error handling here. Type was set to wrong value.
  }

  return val;
}

void BTCSDiffusion::setBoundaryCondition(int index, bctype type, double value) {

  bc[index].type = type;
  bc[index].value = value;
}
