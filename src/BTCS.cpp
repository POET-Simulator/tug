#include <array>
#include <iostream>
#include <tug/BoundaryCondition.hpp>
#include <tug/Diffusion.hpp>
#include <tug/Solver.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/src/Core/Matrix.h>
#include <chrono>
#include <vector>

#include "TugUtils.hpp"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

inline auto
init_delta(const std::array<double, tug::diffusion::MAX_ARR_SIZE> &domain_size,
           const std::array<uint32_t, tug::diffusion::MAX_ARR_SIZE> &grid_cells,
           const uint8_t dim) -> std::vector<double> {
  std::vector<double> out(dim);
  for (uint8_t i = 0; i < dim; i++) {
    // calculate 'size' of each cell in grid
    out[i] = (double)(domain_size.at(i) / grid_cells.at(i)); 
  }
  return out;
}

namespace {
enum { GRID_1D = 1, GRID_2D, GRID_3D };

constexpr int BTCS_MAX_DEP_PER_CELL = 3;
constexpr int BTCS_2D_DT_SIZE = 2;

using DMatrixRowMajor =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using DVectorRowMajor =
    Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>;

inline auto getBCFromFlux(tug::bc::boundary_condition bc, double neighbor_c,
                          double neighbor_alpha) -> double {
  double val = 0;

  if (bc.type == tug::bc::BC_TYPE_CLOSED) {
    val = neighbor_c;
  } else if (bc.type == tug::bc::BC_TYPE_FLUX) {
    // TODO
    //  val = bc[index].value;
  } else {
    // TODO: implement error handling here. Type was set to wrong value.
  }

  return val;
}

auto calc_d_ortho(const DMatrixRowMajor &c, const DMatrixRowMajor &alpha,
                  const tug::bc::BoundaryCondition &bc, bool transposed,
                  double time_step, double dx) -> DMatrixRowMajor {

  uint8_t upper = (transposed ? tug::bc::BC_SIDE_LEFT : tug::bc::BC_SIDE_TOP);
  uint8_t lower =
      (transposed ? tug::bc::BC_SIDE_RIGHT : tug::bc::BC_SIDE_BOTTOM);

  int n_rows = c.rows();
  int n_cols = c.cols();

  DMatrixRowMajor d_ortho(n_rows, n_cols);

  std::array<double, 3> y_values{};

  // first, iterate over first row
  for (int j = 0; j < n_cols; j++) {
    tug::bc::boundary_condition tmp_bc = bc(upper, j);
    double sy = (time_step * alpha(0, j)) / (dx * dx);

    y_values[0] = (tmp_bc.type == tug::bc::BC_TYPE_CONSTANT
                       ? tmp_bc.value
                       : getBCFromFlux(tmp_bc, c(0, j), alpha(0, j)));
    y_values[1] = c(0, j);
    y_values[2] = c(1, j);

    d_ortho(0, j) = -sy * (2 * y_values[0] - 3 * y_values[1] + y_values[2]);
  }

  // then iterate over inlet
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
    tug::bc::boundary_condition tmp_bc = bc(lower, j);
    double sy = (time_step * alpha(end, j)) / (dx * dx);

    y_values[0] = c(end - 1, j);
    y_values[1] = c(end, j);
    y_values[2] = (tmp_bc.type == tug::bc::BC_TYPE_CONSTANT
                       ? tmp_bc.value
                       : getBCFromFlux(tmp_bc, c(end, j), alpha(end, j)));

    d_ortho(end, j) = -sy * (y_values[0] - 3 * y_values[1] + 2 * y_values[2]);
  }

  return d_ortho;
}

auto fillMatrixFromRow(const DVectorRowMajor &alpha,
                       const tug::bc::bc_vec &bc_inner, int size, double dx,
                       double time_step) -> Eigen::SparseMatrix<double> {

  Eigen::SparseMatrix<double> A_matrix(size + 2, size + 2);

  A_matrix.reserve(Eigen::VectorXi::Constant(size + 2, BTCS_MAX_DEP_PER_CELL));

  double sx = 0;

  int A_size = A_matrix.cols();

  A_matrix.insert(0, 0) = 1;

  if (bc_inner[0].type != tug::bc::BC_UNSET) {
    if (bc_inner[0].type != tug::bc::BC_TYPE_CONSTANT) {
      throw_invalid_argument("Inner boundary conditions with other type than "
                             "BC_TYPE_CONSTANT are currently not supported.");
    }
    A_matrix.insert(1, 1) = 1;
  } else {
    sx = (alpha[0] * time_step) / (dx * dx);
    A_matrix.insert(1, 1) = -1. - 3. * sx;
    A_matrix.insert(1, 0) = 2. * sx;
    A_matrix.insert(1, 2) = sx;
  }

  for (int j = 2, k = j - 1; k < size - 1; j++, k++) {
    if (bc_inner[k].type != tug::bc::BC_UNSET) {
      if (bc_inner[k].type != tug::bc::BC_TYPE_CONSTANT) {
        throw_invalid_argument("Inner boundary conditions with other type than "
                               "BC_TYPE_CONSTANT are currently not supported.");
      }
      A_matrix.insert(j, j) = 1;
      continue;
    }
    sx = (alpha[k] * time_step) / (dx * dx);

    A_matrix.insert(j, j) = -1. - 2. * sx;
    A_matrix.insert(j, (j - 1)) = sx;
    A_matrix.insert(j, (j + 1)) = sx;
  }

  if (bc_inner[size - 1].type != tug::bc::BC_UNSET) {
    if (bc_inner[size - 1].type != tug::bc::BC_TYPE_CONSTANT) {
      throw_invalid_argument("Inner boundary conditions with other type than "
                             "BC_TYPE_CONSTANT are currently not supported.");
    }
    A_matrix.insert(A_size - 2, A_size - 2) = 1;
  } else {
    sx = (alpha[size - 1] * time_step) / (dx * dx);
    A_matrix.insert(A_size - 2, A_size - 2) = -1. - 3. * sx;
    A_matrix.insert(A_size - 2, A_size - 3) = sx;
    A_matrix.insert(A_size - 2, A_size - 1) = 2. * sx;
  }

  A_matrix.insert(A_size - 1, A_size - 1) = 1;

  return A_matrix;
}

auto fillVectorFromRow(const DVectorRowMajor &c, const DVectorRowMajor &alpha,
                       const tug::bc::bc_tuple &bc,
                       const tug::bc::bc_vec &bc_inner,
                       const DVectorRowMajor &d_ortho, int size, double dx,
                       double time_step) -> Eigen::VectorXd {

  Eigen::VectorXd b_vector(size + 2);

  tug::bc::boundary_condition left = bc[0];
  tug::bc::boundary_condition right = bc[1];

  bool left_constant = (left.type == tug::bc::BC_TYPE_CONSTANT);
  bool right_constant = (right.type == tug::bc::BC_TYPE_CONSTANT);

  int b_size = b_vector.size();

  for (int j = 0; j < size; j++) {
    if (bc_inner[j].type != tug::bc::BC_UNSET) {
      if (bc_inner[j].type != tug::bc::BC_TYPE_CONSTANT) {
        throw_invalid_argument("Inner boundary conditions with other type than "
                               "BC_TYPE_CONSTANT are currently not supported.");
      }
      b_vector[j + 1] = bc_inner[j].value;
      continue;
    }
    b_vector[j + 1] = -c[j] + d_ortho[j];
  }

  // this is not correct currently.We will fix this when we are able to define
  // FLUX boundary conditions
  b_vector[0] =
      (left_constant ? left.value : getBCFromFlux(left, c[0], alpha[0]));

  b_vector[b_size - 1] =
      (right_constant ? right.value
                      : getBCFromFlux(right, c[size - 1], alpha[size - 1]));

  return b_vector;
}

auto setupBTCSAndSolve(
    DVectorRowMajor &c, const tug::bc::bc_tuple bc_ghost,
    const tug::bc::bc_vec &bc_inner, const DVectorRowMajor &alpha, double dx,
    double time_step, int size, const DVectorRowMajor &d_ortho,
    Eigen::VectorXd (*solver)(const Eigen::SparseMatrix<double> &,
                              const Eigen::VectorXd &)) -> DVectorRowMajor {

  const Eigen::SparseMatrix<double> A_matrix =
      fillMatrixFromRow(alpha, bc_inner, size, dx, time_step);

  const Eigen::VectorXd b_vector = fillVectorFromRow(
      c, alpha, bc_ghost, bc_inner, d_ortho, size, dx, time_step);

  // solving of the LEQ
  Eigen::VectorXd x_vector = solver(A_matrix, b_vector);

  DVectorRowMajor out_vector = x_vector.segment(1, size);

  return out_vector;
}

} // namespace
  //
auto tug::diffusion::BTCS_1D(const tug::diffusion::TugInput &input_param,
                             double *field, const double *alpha) -> double {

  auto start = time_marker();

  uint32_t size = input_param.grid.grid_cells[0];

  auto deltas = init_delta(input_param.grid.domain_size,
                           input_param.grid.grid_cells, GRID_1D);
  double dx = deltas[0];

  double time_step = input_param.time_step;

  const tug::bc::BoundaryCondition bc =
      (input_param.grid.bc != nullptr ? *input_param.grid.bc
                                      : tug::bc::BoundaryCondition(size));

  Eigen::Map<DVectorRowMajor> c_in(field, size);
  Eigen::Map<const DVectorRowMajor> alpha_in(alpha, size);

  DVectorRowMajor input_field = c_in.row(0);

  DVectorRowMajor output = setupBTCSAndSolve(
      input_field, bc.row_boundary(0), bc.getInnerRow(0), alpha_in, dx,
      time_step, size, Eigen::VectorXd::Constant(size, 0), input_param.solver);

  c_in.row(0) << output;

  auto end = time_marker();

  return diff_time(start, end);
}

auto tug::diffusion::ADI_2D(const tug::diffusion::TugInput &input_param,
                            double *field, const double *alpha) -> double {

  auto start = time_marker();

  uint32_t n_cols = input_param.grid.grid_cells[0];
  uint32_t n_rows = input_param.grid.grid_cells[1];

  auto deltas = init_delta(input_param.grid.domain_size,
                           input_param.grid.grid_cells, GRID_2D);
  double dx = deltas[0];
  double dy = deltas[1];

  double local_dt = input_param.time_step / BTCS_2D_DT_SIZE;

  tug::bc::BoundaryCondition bc =
      (input_param.grid.bc != nullptr
           ? *input_param.grid.bc
           : tug::bc::BoundaryCondition(n_cols, n_rows));

  Eigen::Map<DMatrixRowMajor> c_in(field, n_rows, n_cols);
  Eigen::Map<const DMatrixRowMajor> alpha_in(alpha, n_rows, n_cols);

  DMatrixRowMajor d_ortho =
      calc_d_ortho(c_in, alpha_in, bc, false, local_dt, dx);

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n_rows; i++) {
    DVectorRowMajor input_field = c_in.row(i);
    DVectorRowMajor output = setupBTCSAndSolve(
        input_field, bc.row_boundary(i), bc.getInnerRow(i), alpha_in.row(i), dx,
        local_dt, n_cols, d_ortho.row(i), input_param.solver);
    c_in.row(i) << output;
  }

  d_ortho = calc_d_ortho(c_in.transpose(), alpha_in.transpose(), bc, true,
                         local_dt, dy);

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < n_cols; i++) {
    DVectorRowMajor input_field = c_in.col(i);
    DVectorRowMajor output = setupBTCSAndSolve(
        input_field, bc.col_boundary(i), bc.getInnerCol(i), alpha_in.col(i), dy,
        local_dt, n_rows, d_ortho.row(i), input_param.solver);
    c_in.col(i) << output.transpose();
  }

  auto end = time_marker();

  return diff_time(start, end);
}
