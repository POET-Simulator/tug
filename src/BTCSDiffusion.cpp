#include "BTCSDiffusion.hpp"

#include <Eigen/SparseLU>

#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iterator>
#include <ostream>
#include <tuple>
#include <vector>

#include <iostream>

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

void BTCSDiffusion::simulate1D(Eigen::Map<DVectorRowMajor> &c,
                               boundary_condition left,
                               boundary_condition right,
                               const std::vector<double> &alpha, double dx,
                               int size) {

  bool left_is_constant = (left.type == BTCSDiffusion::BC_CONSTANT);
  bool right_is_constant = (right.type == BTCSDiffusion::BC_CONSTANT);

  // The sizes for matrix and vectors of the equation system is defined by the
  // actual size of the input vector and if the system is (partially) closed.
  // Then we will need ghost nodes. So this variable will give the count of
  // ghost nodes.
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

  solveLES();

  // write back result to input/output vector
  c = x_vector.segment(!left_is_constant, c.size());
}

void BTCSDiffusion::simulate2D(Eigen::Map<DMatrixRowMajor> &c,
                               Eigen::Map<const DMatrixRowMajor> &alpha) {

  DMatrixRowMajor tmp_vector;

  int n_cols = c.cols();
  unsigned int size = (this->grid_cells[0] + 2) * (this->grid_cells[1]);

  A_matrix.resize(size, size);
  A_matrix.reserve(Eigen::VectorXi::Constant(size, 3));

  b_vector.resize(size);
  x_vector.resize(size);

  for (int i = 0; i < c.rows(); i++) {
    boundary_condition left = bc[i * n_cols];
    bool left_constant = left.type == BTCSDiffusion::BC_CONSTANT;
    boundary_condition right = bc[((i + 1) * n_cols) - 1];
    bool right_constant = right.type == BTCSDiffusion::BC_CONSTANT;

    fillMatrixFromRow(alpha.row(i), n_cols, i, left_constant, right_constant,
                      deltas[0], this->time_step / 2);
    fillVectorFromRow2D(c, alpha.row(i), i, deltas[0], left, right);
  }

  solveLES();

  tmp_vector = x_vector;
  tmp_vector.transposeInPlace();
  tmp_vector.conservativeResize(c.rows(), c.cols() + 2);

  Eigen::Map<Eigen::MatrixXd> tmp(tmp_vector.data(), c.rows(), c.cols() + 2);

  c = tmp_vector.block(0, 1, c.rows(), c.cols());
  c.transposeInPlace();

  size = (this->grid_cells[0] * (this->grid_cells[1] + 2));

  A_matrix.resize(size, size);
  A_matrix.reserve(Eigen::VectorXi::Constant(size, 3));

  b_vector.resize(size);
  x_vector.resize(size);

  int bottom_offset = bc.size() - (this->grid_cells[0]);
  n_cols = c.cols();

  for (int i = 0; i < c.rows(); i++) {
    boundary_condition left = bc[i];
    bool left_constant = left.type == BTCSDiffusion::BC_CONSTANT;
    boundary_condition right = bc[bottom_offset + i];
    bool right_constant = right.type == BTCSDiffusion::BC_CONSTANT;

    fillMatrixFromRow(alpha.col(i), n_cols, i, left_constant, right_constant,
                      deltas[1], this->time_step / 2);
    fillVectorFromRow2D(c, alpha.row(i), i, deltas[1], left, right);
  }

  solveLES();

  tmp_vector = x_vector;
  tmp_vector.transposeInPlace();
  tmp_vector.conservativeResize(c.rows(), c.cols() + 2);

  c = tmp_vector.block(0, 1, c.rows(), c.cols());

  c.transposeInPlace();
}

void BTCSDiffusion::fillMatrixFromRow(const DVectorRowMajor &alpha, int n_cols,
                                      int row, bool left_constant,
                                      bool right_constant, double delta,
                                      double time_step) {

  n_cols += 2;
  int offset = n_cols * row;

  A_matrix.insert(offset, offset) = 1;

  if (left_constant)
    A_matrix.insert(offset + 1, offset + 1) = 1;

  A_matrix.insert(offset + (n_cols - 1), offset + (n_cols - 1)) = 1;

  if (right_constant)
    A_matrix.insert(offset + (n_cols - 2), offset + (n_cols - 2)) = 1;

  for (int j = 1 + left_constant, k = j - 1; j < n_cols - (1 - right_constant);
       j++, k++) {
    double sx = (alpha[j - 1] * time_step) / (delta * delta);

    if (this->bc[row * (n_cols - 2) + k].type == BTCSDiffusion::BC_CONSTANT) {
      A_matrix.insert(offset + j, offset + j) = 1;
      continue;
    }

    A_matrix.insert(offset + j, offset + j) = -1. - 2. * sx;
    A_matrix.insert(offset + j, offset + (j - 1)) = sx;
    A_matrix.insert(offset + j, offset + (j + 1)) = sx;
  }
}

void BTCSDiffusion::fillVectorFromRow2D(Eigen::Map<DMatrixRowMajor> &c,
                                        const Eigen::VectorXd alpha, int row,
                                        double delta, boundary_condition left,
                                        boundary_condition right) {

  int ncol = c.cols();
  int nrow = c.rows();
  int offset = ncol + 2;

  if (left.type != BTCSDiffusion::BC_CONSTANT) {
    // this is not correct currently.We will fix this when we are able to define
    // FLUX boundary conditions
    b_vector[offset * row] = getBCFromFlux(left, c(row, 0), alpha[0]);
  }

  if (right.type != BTCSDiffusion::BC_CONSTANT) {
    b_vector[offset * row + (offset - 1)] =
        getBCFromFlux(right, c(row, ncol - 1), alpha[ncol - 1]);
  }

  for (int j = 1; j < offset - 1; j++) {
    boundary_condition tmp_bc = this->bc[ncol * row + (j - 1)];

    if (tmp_bc.type == BTCSDiffusion::BC_CONSTANT) {
      b_vector[offset * row + j] = tmp_bc.value;
    } else {

      double y_values[3];
      y_values[0] =
          (row != 0 ? c(row - 1, j - 1)
                    : getBCFromFlux(tmp_bc, c(row, j - 1), alpha[j - 1]));
      y_values[1] = c(row, j - 1);
      y_values[2] = (row != nrow - 1
                         ? c(row + 1, j - 1)
                         : getBCFromFlux(tmp_bc, c(row, j - 1), alpha[j - 1]));

      double t0_c =
          alpha[j - 1] *
          ((y_values[0] - 2 * y_values[1] + y_values[2]) / (delta * delta));
      b_vector[offset * row + j] = -c(row, j - 1) - t0_c;
    }
  }
}

void BTCSDiffusion::setTimestep(double time_step) {
  this->time_step = time_step;
}

void BTCSDiffusion::simulate(std::vector<double> &c,
                             const std::vector<double> &alpha) {
  if (this->grid_dim == 1) {
    assert(c.size() == grid_cells[0]);
    Eigen::Map<DVectorRowMajor> c_in(c.data(), this->grid_cells[0]);
    simulate1D(c_in, bc[0], bc[grid_cells[0] + 1], alpha, this->deltas[0],
               this->grid_cells[0]);
  }
  if (this->grid_dim == 2) {
    assert(c.size() == grid_cells[0] * grid_cells[1]);
    Eigen::Map<DMatrixRowMajor> c_in(c.data(), this->grid_cells[1],
                                     this->grid_cells[0]);

    Eigen::Map<const DMatrixRowMajor> alpha_in(
        alpha.data(), this->grid_cells[1], this->grid_cells[0]);

    simulate2D(c_in, alpha_in);
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

inline void BTCSDiffusion::solveLES() {
  // start to solve
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A_matrix);

  solver.factorize(A_matrix);

  x_vector = solver.solve(b_vector);
}
