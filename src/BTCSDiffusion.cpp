#include "diffusion/BTCSDiffusion.hpp"
#include "diffusion/BTCSBoundaryCondition.hpp"

#include <Eigen/SparseLU>

#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <Eigen/src/SparseCore/SparseMatrixBase.h>
#include <algorithm>
#include <array>
#include <bits/stdint-uintn.h>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <ostream>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

constexpr double DOUBLE_MACHINE_EPSILON = 1.93e-34;

constexpr int BTCS_MAX_DEP_PER_CELL = 3;
constexpr int BTCS_2D_DT_SIZE = 2;
Diffusion::BTCSDiffusion::BTCSDiffusion(unsigned int dim) : grid_dim(dim) {

  grid_cells.resize(dim, 1);
  domain_size.resize(dim, 1);
  deltas.resize(dim, 1);

  this->time_step = 0;
}

void Diffusion::BTCSDiffusion::setXDimensions(double domain_size,
                                              unsigned int n_grid_cells) {
  this->domain_size[0] = domain_size;
  this->grid_cells[0] = n_grid_cells;

  updateInternals();
}

void Diffusion::BTCSDiffusion::setYDimensions(double domain_size,
                                              unsigned int n_grid_cells) {
  this->domain_size[1] = domain_size;
  this->grid_cells[1] = n_grid_cells;

  updateInternals();
}

void Diffusion::BTCSDiffusion::setZDimensions(double domain_size,
                                              unsigned int n_grid_cells) {
  this->domain_size[2] = domain_size;
  this->grid_cells[2] = n_grid_cells;

  updateInternals();
}

auto Diffusion::BTCSDiffusion::getXGridCellsN() -> unsigned int {
  return this->grid_cells[0];
}
auto Diffusion::BTCSDiffusion::getYGridCellsN() -> unsigned int {
  return this->grid_cells[1];
}
auto Diffusion::BTCSDiffusion::getZGridCellsN() -> unsigned int {
  return this->grid_cells[2];
}
auto Diffusion::BTCSDiffusion::getXDomainSize() -> double {
  return this->domain_size[0];
}
auto Diffusion::BTCSDiffusion::getYDomainSize() -> double {
  return this->domain_size[1];
}
auto Diffusion::BTCSDiffusion::getZDomainSize() -> double {
  return this->domain_size[2];
}

void Diffusion::BTCSDiffusion::updateInternals() {
  for (int i = 0; i < grid_dim; i++) {
    deltas[i] = (double)domain_size[i] / grid_cells[i];
  }
}
void Diffusion::BTCSDiffusion::simulate_base(DVectorRowMajor &c,
                                             const Diffusion::bc_tuple &bc,
                                             const DVectorRowMajor &alpha,
                                             double dx, double time_step,
                                             int size,
                                             const DVectorRowMajor &d_ortho) {

  Eigen::SparseMatrix<double> A_matrix;
  Eigen::VectorXd b_vector;
  Eigen::VectorXd x_vector;

  A_matrix.resize(size + 2, size + 2);
  A_matrix.reserve(Eigen::VectorXi::Constant(size + 2, BTCS_MAX_DEP_PER_CELL));

  b_vector.resize(size + 2);
  x_vector.resize(size + 2);

  fillMatrixFromRow(A_matrix, alpha, size, dx, time_step);
  fillVectorFromRow(b_vector, c, alpha, bc, d_ortho, size, dx, time_step);

  // start to solve
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(A_matrix);

  solver.factorize(A_matrix);

  x_vector = solver.solve(b_vector);

  c = x_vector.segment(1, size);
}

void Diffusion::BTCSDiffusion::simulate1D(
    Eigen::Map<DVectorRowMajor> &c, Eigen::Map<const DVectorRowMajor> &alpha,
    const Diffusion::BTCSBoundaryCondition &bc) {

  int size = this->grid_cells[0];
  double dx = this->deltas[0];
  double time_step = this->time_step;

  DVectorRowMajor input_field = c.row(0);

  simulate_base(input_field, bc.row(0), alpha, dx, time_step, size,
                Eigen::VectorXd::Constant(size, 0));

  c.row(0) << input_field;
}

void Diffusion::BTCSDiffusion::simulate2D(
    Eigen::Map<DMatrixRowMajor> &c, Eigen::Map<const DMatrixRowMajor> &alpha,
    const Diffusion::BTCSBoundaryCondition &bc) {

  int n_rows = this->grid_cells[1];
  int n_cols = this->grid_cells[0];
  double dx = this->deltas[0];

  double local_dt = this->time_step / BTCS_2D_DT_SIZE;

  DMatrixRowMajor d_ortho = calc_d_ortho(c, alpha, bc, false, local_dt, dx);

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n_rows; i++) {
    DVectorRowMajor input_field = c.row(i);
    simulate_base(input_field, bc.row(i), alpha.row(i), dx, local_dt, n_cols,
                  d_ortho.row(i));
    c.row(i) << input_field;
  }

  dx = this->deltas[1];

  d_ortho =
      calc_d_ortho(c.transpose(), alpha.transpose(), bc, true, local_dt, dx);

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n_cols; i++) {
    DVectorRowMajor input_field = c.col(i);
    simulate_base(input_field, bc.col(i), alpha.col(i), dx, local_dt, n_rows,
                  d_ortho.row(i));
    c.col(i) << input_field.transpose();
  }
}

auto Diffusion::BTCSDiffusion::calc_d_ortho(
    const DMatrixRowMajor &c, const DMatrixRowMajor &alpha,
    const Diffusion::BTCSBoundaryCondition &bc, bool transposed,
    double time_step, double dx) -> DMatrixRowMajor {

  uint8_t upper = (transposed ? Diffusion::BC_SIDE_LEFT : Diffusion::BC_SIDE_TOP);
  uint8_t lower = (transposed ? Diffusion::BC_SIDE_RIGHT : Diffusion::BC_SIDE_BOTTOM);

  int n_rows = c.rows();
  int n_cols = c.cols();

  DMatrixRowMajor d_ortho(n_rows, n_cols);

  std::array<double, 3> y_values{};

  // first, iterate over first row
  for (int j = 0; j < n_cols; j++) {
    boundary_condition tmp_bc = bc(upper, j);
    double sy = (time_step * alpha(0, j)) / (dx * dx);

    y_values[0] = (tmp_bc.type == Diffusion::BC_TYPE_CONSTANT
                       ? tmp_bc.value
                       : getBCFromFlux(tmp_bc, c(0, j), alpha(0, j)));
    y_values[1] = c(0, j);
    y_values[2] = c(1, j);

    d_ortho(0, j) = -sy * (2 * y_values[0] - 3 * y_values[1] + y_values[2]);
  }

// then iterate over inlet
#pragma omp parallel for private(y_values) schedule(dynamic)
  for (int i = 1; i < n_rows - 1; i++) {
    for (int j = 0; j < n_cols; j++) {
      double sy = (time_step * alpha(i, j)) / (dx * dx);

      y_values[0] = c(i - 1, j);
      y_values[1] = c(i, j);
      y_values[2] = c(i + 1, j);

      d_ortho(i, j) = -sy * (y_values[0] - 2 * y_values[1] + y_values[2]);
    }
  }

  int end = n_rows - 1;

  // and finally over last row
  for (int j = 0; j < n_cols; j++) {
    boundary_condition tmp_bc = bc(lower, j);
    double sy = (time_step * alpha(end, j)) / (dx * dx);

    y_values[0] = c(end - 1, j);
    y_values[1] = c(end, j);
    y_values[2] = (tmp_bc.type == Diffusion::BC_TYPE_CONSTANT
                       ? tmp_bc.value
                       : getBCFromFlux(tmp_bc, c(end, j), alpha(end, j)));

    d_ortho(end, j) = -sy * (y_values[0] - 3 * y_values[1] + 2 * y_values[2]);
  }

  return d_ortho;
}

void Diffusion::BTCSDiffusion::fillMatrixFromRow(
    Eigen::SparseMatrix<double> &A_matrix, const DVectorRowMajor &alpha,
    int size, double dx, double time_step) {

  double sx = 0;

  int A_size = A_matrix.cols();

  A_matrix.insert(0, 0) = 1;

  sx = (alpha[0] * time_step) / (dx * dx);
  A_matrix.insert(1, 1) = -1. - 3. * sx;
  A_matrix.insert(1, 0) = 2. * sx;
  A_matrix.insert(1, 2) = sx;

  for (int j = 2, k = j - 1; k < size - 1; j++, k++) {
    sx = (alpha[k] * time_step) / (dx * dx);

    A_matrix.insert(j, j) = -1. - 2. * sx;
    A_matrix.insert(j, (j - 1)) = sx;
    A_matrix.insert(j, (j + 1)) = sx;
  }

  sx = (alpha[size - 1] * time_step) / (dx * dx);
  A_matrix.insert(A_size - 2, A_size - 2) = -1. - 3. * sx;
  A_matrix.insert(A_size - 2, A_size - 3) = sx;
  A_matrix.insert(A_size - 2, A_size - 1) = 2. * sx;

  A_matrix.insert(A_size - 1, A_size - 1) = 1;
}

void Diffusion::BTCSDiffusion::fillVectorFromRow(
    Eigen::VectorXd &b_vector, const DVectorRowMajor &c,
    const DVectorRowMajor &alpha, const bc_tuple &bc,
    const DVectorRowMajor &d_ortho, int size, double dx, double time_step) {

  Diffusion::boundary_condition left = bc[0];
  Diffusion::boundary_condition right = bc[1];

  bool left_constant = (left.type == Diffusion::BC_TYPE_CONSTANT);
  bool right_constant = (right.type == Diffusion::BC_TYPE_CONSTANT);

  int b_size = b_vector.size();

  for (int j = 0; j < size; j++) {
    b_vector[j + 1] = -c[j] + d_ortho[j];
  }

  // this is not correct currently.We will fix this when we are able to define
  // FLUX boundary conditions
  b_vector[0] =
      (left_constant ? left.value : getBCFromFlux(left, c[0], alpha[0]));

  b_vector[b_size - 1] =
      (right_constant ? right.value
                      : getBCFromFlux(right, c[size - 1], alpha[size - 1]));
}

void Diffusion::BTCSDiffusion::setTimestep(double time_step) {
  this->time_step = time_step;
}

auto Diffusion::BTCSDiffusion::simulate(
    double *c, double *alpha, const Diffusion::BTCSBoundaryCondition &bc)
    -> double {

  std::chrono::high_resolution_clock::time_point start =
      std::chrono::high_resolution_clock::now();

  if (this->grid_dim == 1) {
    Eigen::Map<DVectorRowMajor> c_in(c, this->grid_cells[0]);
    Eigen::Map<const DVectorRowMajor> alpha_in(alpha, this->grid_cells[0]);

    simulate1D(c_in, alpha_in, bc);
  }
  if (this->grid_dim == 2) {
    Eigen::Map<DMatrixRowMajor> c_in(c, this->grid_cells[1],
                                     this->grid_cells[0]);

    Eigen::Map<const DMatrixRowMajor> alpha_in(alpha, this->grid_cells[1],
                                               this->grid_cells[0]);

    simulate2D(c_in, alpha_in, bc);
  }

  std::chrono::high_resolution_clock::time_point end =
      std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> duration =
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

  return duration.count();
}

inline auto Diffusion::BTCSDiffusion::getBCFromFlux(boundary_condition bc,
                                                    double neighbor_c,
                                                    double neighbor_alpha)
    -> double {

  double val = 0;

  if (bc.type == Diffusion::BC_TYPE_CLOSED) {
    val = neighbor_c;
  } else if (bc.type == Diffusion::BC_TYPE_FLUX) {
    // TODO
    //  val = bc[index].value;
  } else {
    // TODO: implement error handling here. Type was set to wrong value.
  }

  return val;
}
