#include "BTCSDiffusion.hpp"
#include "BoundaryCondition.hpp"

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

#define BTCS_MAX_DEP_PER_CELL 3

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

void Diffusion::BTCSDiffusion::simulate_base(
    DVectorRowMajor &c, Eigen::Map<const BCVectorRowMajor> &bc,
    Eigen::Map<const DVectorRowMajor> &alpha, double dx, double time_step,
    int size, DVectorRowMajor &t0_c) {

  // The sizes for matrix and vectors of the equation system is defined by the
  // actual size of the input vector and if the system is (partially) closed.
  // Then we will need ghost nodes. So this variable will give the count of
  // ghost nodes.
  // int bc_offset = !left_is_constant + !right_is_constant;
  // ;

  // set sizes of private and yet allocated vectors
  // b_vector.resize(size + bc_offset);
  // x_vector.resize(size + bc_offset);

  // /*
  //  * Begin to solve the equation system using LU solver of Eigen.
  //  *
  //  * But first fill the A matrix and b vector.
  //  */

  // // Set boundary condition for ghost nodes (for closed or flux system) or
  // outer
  // // inlet nodes (constant boundary condition)
  // A_matrix.resize(size + bc_offset, size + bc_offset);
  // A_matrix.reserve(Eigen::VectorXi::Constant(size + bc_offset, 3));

  // A_matrix.insert(0, 0) = 1;
  // b_vector[0] =
  //     (left_is_constant ? left.value : getBCFromFlux(left, c[0], alpha[0]));

  // A_matrix.insert(size + 1, size + 1) = 1;
  // b_vector[size + 1] =
  //     (right_is_constant ? right.value
  //                        : getBCFromFlux(right, c[size - 1], alpha[size -
  //                        1]));

  // Start filling the A matrix
  // =i= is used for equation system matrix and vector indexing
  // and =j= for indexing of c,alpha and bc
  // for (int i = 1, j = i + !(left_is_constant); i < size - right_is_constant;
  //      i++, j++) {

  //   // if current grid cell is considered as constant boundary conditon
  //   if (bc[j].type == Diffusion::BC_CONSTANT) {
  //     A_matrix.insert(i, i) = 1;
  //     b_vector[i] = bc[j].value;
  //     continue;
  //   }

  //   double sx = (alpha[j] * time_step) / (dx * dx);

  //   A_matrix.insert(i, i) = -1. - 2. * sx;
  //   A_matrix.insert(i, i - 1) = sx;
  //   A_matrix.insert(i, i + 1) = sx;

  //   b_vector[i] = -c[j];
  // }

  fillMatrixFromRow(alpha, bc, size, dx, time_step);
  fillVectorFromRowADI(c, alpha, bc, t0_c, size, dx, time_step);

  solveLES();

  // write back result to input/output vector
  // c = x_vector.segment(!left_is_constant, c.size());
}

inline void Diffusion::BTCSDiffusion::reserveMemory(int size,
                                                    int max_count_per_line) {
  size += 2;

  A_matrix.resize(size, size);
  A_matrix.reserve(Eigen::VectorXi::Constant(size, max_count_per_line));

  b_vector.resize(size);
  x_vector.resize(size);
}

void Diffusion::BTCSDiffusion::simulate1D(
    Eigen::Map<DVectorRowMajor> &c, Eigen::Map<const DVectorRowMajor> &alpha,
    Eigen::Map<const BCVectorRowMajor> &bc) {

  int size = this->grid_cells[0];
  double dx = this->deltas[0];
  double time_step = this->time_step;

  reserveMemory(size, BTCS_MAX_DEP_PER_CELL);

  fillMatrixFromRow(alpha.row(0), bc.row(0), size, dx, time_step);
  fillVectorFromRowADI(c, alpha, bc, Eigen::VectorXd::Constant(size, 0), size,
                       dx, time_step);

  std::cout << A_matrix << std::endl;

  solveLES();

  c = x_vector.segment(1, size);
}

void Diffusion::BTCSDiffusion::simulate2D(
    Eigen::Map<DMatrixRowMajor> &c, Eigen::Map<const DMatrixRowMajor> &alpha,
    Eigen::Map<const BCMatrixRowMajor> &bc) {

  // double local_dt = this->time_step / 2.;
  // DMatrixRowMajor tmp_vector;

  // int n_cols = c.cols();
  // unsigned int size = (this->grid_cells[0] + 2) * (this->grid_cells[1]);

  // A_matrix.resize(size, size);
  // A_matrix.reserve(Eigen::VectorXi::Constant(size, 3));

  // b_vector.resize(size);
  // x_vector.resize(size);

  // for (int i = 0; i < c.rows(); i++) {
  //   boundary_condition left = bc(i, 0);
  //   bool left_constant = left.type == Diffusion::BC_CONSTANT;
  //   boundary_condition right = bc(i, n_cols - 1);
  //   bool right_constant = right.type == Diffusion::BC_CONSTANT;

  //   fillMatrixFromRow(alpha.row(i), n_cols, i, left_constant, right_constant,
  //                     deltas[0], this->time_step / 2, bc.row(i));
  //   fillVectorFromRowADI(c, alpha.row(i), i, deltas[0], left, right,
  //   local_dt,
  //                        bc.row(i));
  // }

  // solveLES();

  // tmp_vector = x_vector;
  // tmp_vector.transposeInPlace();
  // tmp_vector.conservativeResize(c.rows(), c.cols() + 2);

  // Eigen::Map<Eigen::MatrixXd> tmp(tmp_vector.data(), c.rows(), c.cols() + 2);

  // c = tmp_vector.block(0, 1, c.rows(), c.cols());
  // c.transposeInPlace();

  // size = (this->grid_cells[0] * (this->grid_cells[1] + 2));

  // A_matrix.resize(size, size);
  // A_matrix.reserve(Eigen::VectorXi::Constant(size, 3));

  // b_vector.resize(size);
  // x_vector.resize(size);

  // n_cols = c.cols();

  // for (int i = 0; i < c.rows(); i++) {
  //   boundary_condition left = bc(0, i);
  //   bool left_constant = left.type == Diffusion::BC_CONSTANT;
  //   boundary_condition right = bc(n_cols - 1, i);
  //   bool right_constant = right.type == Diffusion::BC_CONSTANT;

  //   fillMatrixFromRow(alpha.col(i), n_cols, i, left_constant, right_constant,
  //                     deltas[1], this->time_step / 2, bc.col(i));
  //   fillVectorFromRowADI(c, alpha.row(i), i, deltas[1], left, right,
  //   local_dt,
  //                        bc.col(i));
  // }

  // solveLES();

  // tmp_vector = x_vector;
  // tmp_vector.transposeInPlace();
  // tmp_vector.conservativeResize(c.rows(), c.cols() + 2);

  // c = tmp_vector.block(0, 1, c.rows(), c.cols());

  // c.transposeInPlace();
}

inline void Diffusion::BTCSDiffusion::fillMatrixFromRow(
    const Eigen::VectorXd &alpha,
    const Eigen::Vector<Diffusion::boundary_condition, Eigen::Dynamic> &bc,
    int size, double dx, double time_step) {

  Diffusion::boundary_condition left = bc[0];
  Diffusion::boundary_condition right = bc[size - 1];

  bool left_constant = (left.type == Diffusion::BC_CONSTANT);
  bool right_constant = (right.type == Diffusion::BC_CONSTANT);

  int A_size = A_matrix.cols();

  A_matrix.insert(0, 0) = 1;

  if (left_constant)
    A_matrix.insert(1, 1) = 1;

  A_matrix.insert(A_size - 1, A_size - 1) = 1;

  if (right_constant)
    A_matrix.insert(A_size - 2, A_size - 2) = 1;

  for (int j = 1 + left_constant, k = j - 1; j < size - (1 - right_constant);
       j++, k++) {
    double sx = (alpha[k] * time_step) / (dx * dx);

    if (bc[k].type == Diffusion::BC_CONSTANT) {
      A_matrix.insert(j, j) = 1;
      continue;
    }

    A_matrix.insert(j, j) = -1. - 2. * sx;
    A_matrix.insert(j, (j - 1)) = sx;
    A_matrix.insert(j, (j + 1)) = sx;
  }
}

inline void Diffusion::BTCSDiffusion::fillVectorFromRowADI(
    const DVectorRowMajor &c, const Eigen::VectorXd alpha,
    const BCVectorRowMajor &bc, const DVectorRowMajor &t0_c, int size,
    double dx, double time_step) {

  Diffusion::boundary_condition left = bc[0];
  Diffusion::boundary_condition right = bc[size - 1];

  bool left_constant = (left.type == Diffusion::BC_CONSTANT);
  bool right_constant = (right.type == Diffusion::BC_CONSTANT);

  int b_size = b_vector.size();

  for (int j = 0; j < size; j++) {
    boundary_condition tmp_bc = bc[j];

    if (tmp_bc.type == Diffusion::BC_CONSTANT) {
      b_vector[j + 1] = tmp_bc.value;
      continue;
    }

    // double y_values[3];
    // y_values[0] =
    //     (row != 0 ? c(row - 1, j) : getBCFromFlux(tmp_bc, c(row, j),
    //     alpha[j]));
    // y_values[1] = c(row, j);
    // y_values[2] =
    //     (row != nrow - 1 ? c(row + 1, j)
    //                      : getBCFromFlux(tmp_bc, c(row, j), alpha[j]));

    double t0_c_j = time_step * alpha[j] * (t0_c[j] / (dx * dx));
    b_vector[j + 1] = -c[j] - t0_c_j;
  }

  if (!left_constant) {
    // this is not correct currently.We will fix this when we are able to define
    // FLUX boundary conditions
    b_vector[0] = getBCFromFlux(left, b_vector[1], alpha[0]);
  }

  if (!right_constant) {
    b_vector[b_size - 1] =
        getBCFromFlux(right, b_vector[size - 2], alpha[size - 1]);
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

    simulate1D(c_in, alpha_in, bc_in);
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
