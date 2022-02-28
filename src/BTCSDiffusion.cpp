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

Diffusion::BTCSDiffusion::BTCSDiffusion(unsigned int dim) : grid_dim(dim) {
  assert(dim <= 3);

  grid_cells.resize(dim, 1);
  domain_size.resize(dim, 1);
  deltas.resize(dim, 1);
}

void Diffusion::BTCSDiffusion::setXDimensions(double domain_size,
                                              unsigned int n_grid_cells) {
  assert(this->grid_dim > 0);
  this->domain_size[0] = domain_size;
  this->grid_cells[0] = n_grid_cells;

  updateInternals();
}

void Diffusion::BTCSDiffusion::setYDimensions(double domain_size,
                                              unsigned int n_grid_cells) {
  assert(this->grid_dim > 1);
  this->domain_size[1] = domain_size;
  this->grid_cells[1] = n_grid_cells;

  updateInternals();
}

void Diffusion::BTCSDiffusion::setZDimensions(double domain_size,
                                              unsigned int n_grid_cells) {
  assert(this->grid_dim > 2);
  this->domain_size[2] = domain_size;
  this->grid_cells[2] = n_grid_cells;

  updateInternals();
}

unsigned int Diffusion::BTCSDiffusion::getXGridCellsN() {
  return this->grid_cells[0];
}
unsigned int Diffusion::BTCSDiffusion::getYGridCellsN() {
  return this->grid_cells[1];
}
unsigned int Diffusion::BTCSDiffusion::getZGridCellsN() {
  return this->grid_cells[2];
}
unsigned int Diffusion::BTCSDiffusion::getXDomainSize() {
  return this->domain_size[0];
}
unsigned int Diffusion::BTCSDiffusion::getYDomainSize() {
  return this->domain_size[1];
}
unsigned int Diffusion::BTCSDiffusion::getZDomainSize() {
  return this->domain_size[2];
}

void Diffusion::BTCSDiffusion::updateInternals() {
  for (int i = 0; i < grid_dim; i++) {
    deltas[i] = (double)domain_size[i] / grid_cells[i];
  }
}

void Diffusion::BTCSDiffusion::simulate1D(
    Eigen::Map<DVectorRowMajor> &c, Diffusion::boundary_condition left,
    Diffusion::boundary_condition right, Eigen::Map<const BCVectorRowMajor> &bc,
    Eigen::Map<const DVectorRowMajor> &alpha, double dx, int size) {

  bool left_is_constant = (left.type == Diffusion::BC_CONSTANT);
  bool right_is_constant = (right.type == Diffusion::BC_CONSTANT);

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
    if (bc[j].type == Diffusion::BC_CONSTANT) {
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

void Diffusion::BTCSDiffusion::simulate2D(
    Eigen::Map<DMatrixRowMajor> &c, Eigen::Map<const DMatrixRowMajor> &alpha,
    Eigen::Map<const BCMatrixRowMajor> &bc) {

  double local_dt = this->time_step / 2.;
  DMatrixRowMajor tmp_vector;

  int n_cols = c.cols();
  unsigned int size = (this->grid_cells[0] + 2) * (this->grid_cells[1]);

  A_matrix.resize(size, size);
  A_matrix.reserve(Eigen::VectorXi::Constant(size, 3));

  b_vector.resize(size);
  x_vector.resize(size);

  for (int i = 0; i < c.rows(); i++) {
    boundary_condition left = bc(i, 0);
    bool left_constant = left.type == Diffusion::BC_CONSTANT;
    boundary_condition right = bc(i, n_cols - 1);
    bool right_constant = right.type == Diffusion::BC_CONSTANT;

    fillMatrixFromRow(alpha.row(i), n_cols, i, left_constant, right_constant,
                      deltas[0], this->time_step / 2, bc.row(i));
    fillVectorFromRowADI(c, alpha.row(i), i, deltas[0], left, right, local_dt, bc.row(i));
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

  n_cols = c.cols();

  for (int i = 0; i < c.rows(); i++) {
    boundary_condition left = bc(0, i);
    bool left_constant = left.type == Diffusion::BC_CONSTANT;
    boundary_condition right = bc(n_cols - 1, i);
    bool right_constant = right.type == Diffusion::BC_CONSTANT;

    fillMatrixFromRow(alpha.col(i), n_cols, i, left_constant, right_constant,
                      deltas[1], this->time_step / 2, bc.col(i));
    fillVectorFromRowADI(c, alpha.row(i), i, deltas[1], left, right, local_dt, bc.col(i));
  }

  solveLES();

  tmp_vector = x_vector;
  tmp_vector.transposeInPlace();
  tmp_vector.conservativeResize(c.rows(), c.cols() + 2);

  c = tmp_vector.block(0, 1, c.rows(), c.cols());

  c.transposeInPlace();
}

void Diffusion::BTCSDiffusion::fillMatrixFromRow(const DVectorRowMajor &alpha,
                                                 int n_cols, int row,
                                                 bool left_constant,
                                                 bool right_constant,
                                                 double delta, double time_step,
                                                 const BCVectorRowMajor &bc) {

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

    if (bc[k].type == Diffusion::BC_CONSTANT) {
      A_matrix.insert(offset + j, offset + j) = 1;
      continue;
    }

    A_matrix.insert(offset + j, offset + j) = -1. - 2. * sx;
    A_matrix.insert(offset + j, offset + (j - 1)) = sx;
    A_matrix.insert(offset + j, offset + (j + 1)) = sx;
  }
}

void Diffusion::BTCSDiffusion::fillVectorFromRowADI(
    Eigen::Map<DMatrixRowMajor> &c, const Eigen::VectorXd alpha, int row,
    double delta, boundary_condition left, boundary_condition right,
    double time_step, const BCVectorRowMajor &bc) {

  int ncol = c.cols();
  int nrow = c.rows();
  int offset = ncol + 2;

  if (left.type != Diffusion::BC_CONSTANT) {
    // this is not correct currently.We will fix this when we are able to define
    // FLUX boundary conditions
    b_vector[offset * row] = getBCFromFlux(left, c(row, 0), alpha[0]);
  }

  if (right.type != Diffusion::BC_CONSTANT) {
    b_vector[offset * row + (offset - 1)] =
        getBCFromFlux(right, c(row, ncol - 1), alpha[ncol - 1]);
  }

  for (int j = 0; j < ncol; j++) {
    boundary_condition tmp_bc = bc[j];

    if (tmp_bc.type == Diffusion::BC_CONSTANT) {
      b_vector[offset * row + (j + 1)] = tmp_bc.value;
      continue;
    }

    double y_values[3];
    y_values[0] =
        (row != 0 ? c(row - 1, j) : getBCFromFlux(tmp_bc, c(row, j), alpha[j]));
    y_values[1] = c(row, j);
    y_values[2] =
        (row != nrow - 1 ? c(row + 1, j)
                         : getBCFromFlux(tmp_bc, c(row, j), alpha[j]));

    double t0_c =
        time_step * alpha[j] *
        ((y_values[0] - 2 * y_values[1] + y_values[2]) / (delta * delta));
    b_vector[offset * row + (j + 1)] = -c(row, j) - (t0_c);
  }
}

void Diffusion::BTCSDiffusion::setTimestep(double time_step) {
  this->time_step = time_step;
}

void Diffusion::BTCSDiffusion::simulate(double *c, double *alpha,
                                        Diffusion::boundary_condition *bc) {
  if (this->grid_dim == 1) {
    Eigen::Map<DVectorRowMajor> c_in(c, this->grid_cells[0]);
    Eigen::Map<const DVectorRowMajor> alpha_in(alpha, this->grid_cells[0]);
    Eigen::Map<const BCVectorRowMajor> bc_in(bc, this->grid_cells[0]);

    simulate1D(c_in, bc[0], bc[this->grid_cells[0] - 1], bc_in, alpha_in,
               this->deltas[0], this->grid_cells[0]);
  }
  if (this->grid_dim == 2) {
    Eigen::Map<DMatrixRowMajor> c_in(c, this->grid_cells[1],
                                     this->grid_cells[0]);

    Eigen::Map<const DMatrixRowMajor> alpha_in(alpha, this->grid_cells[1],
                                               this->grid_cells[0]);

    Eigen::Map<const BCMatrixRowMajor> bc_in(bc, this->grid_cells[1],
                                             this->grid_cells[0]);

    simulate2D(c_in, alpha_in, bc_in);
  }
}

inline double Diffusion::BTCSDiffusion::getBCFromFlux(boundary_condition bc,
                                                      double neighbor_c,
                                                      double neighbor_alpha) {

  double val;

  if (bc.type == Diffusion::BC_CLOSED) {
    val = neighbor_c;
  } else if (bc.type == Diffusion::BC_FLUX) {
    // TODO
    //  val = bc[index].value;
  } else {
    // TODO: implement error handling here. Type was set to wrong value.
  }

  return val;
}

inline void Diffusion::BTCSDiffusion::solveLES() {
  // start to solve
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A_matrix);

  solver.factorize(A_matrix);

  x_vector = solver.solve(b_vector);
}
