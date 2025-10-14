/**
 * @file BTCS.hpp
 * @brief Implementation of heterogenous BTCS (backward time-centered space)
 * solution of diffusion equation in 1D and 2D space. Internally the
 * alternating-direction implicit (ADI) method is used. Version 2, because
 * Version 1 was an implementation for the homogeneous BTCS solution.
 *
 */

#ifndef BTCS_H_
#define BTCS_H_

#include <cstddef>
#include <tug/Boundary.hpp>
#include <tug/Core/Matrix.hpp>
#include <tug/Core/Numeric/SimulationInput.hpp>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

namespace tug {

// optimization to remove Eigen sparse matrix
template <class T> class Diagonals {
public:
  Diagonals() : left(), center(), right() {};
  Diagonals(std::size_t size) : left(size), center(size), right(size) {};

public:
  std::vector<T> left;
  std::vector<T> center;
  std::vector<T> right;
};

// calculates coefficient for boundary in constant case
template <class T>
constexpr std::pair<T, T> calcBoundaryCoeffConstant(T alpha_center,
                                                    T alpha_side, T sx) {
  const T centerCoeff = 1 + sx * (calcAlphaIntercell(alpha_center, alpha_side) +
                                  2 * alpha_center);
  const T sideCoeff = -sx * calcAlphaIntercell(alpha_center, alpha_side);

  return {centerCoeff, sideCoeff};
}

// calculates coefficient for boundary in closed case
template <class T>
constexpr std::pair<T, T> calcBoundaryCoeffClosed(T alpha_center, T alpha_side,
                                                  T sx) {
  const T centerCoeff = 1 + sx * calcAlphaIntercell(alpha_center, alpha_side);

  const T sideCoeff = -sx * calcAlphaIntercell(alpha_center, alpha_side);

  return {centerCoeff, sideCoeff};
}

// creates coefficient matrix for next time step from alphas in x-direction
template <class T>
static Diagonals<T>
createCoeffMatrix(const RowMajMat<T> &alpha,
                  const std::vector<BoundaryElement<T>> &bcLeft,
                  const std::vector<BoundaryElement<T>> &bcRight,
                  const std::vector<std::pair<bool, T>> &inner_bc, int numCols,
                  int rowIndex, T sx) {

  // square matrix of column^2 dimension for the coefficients
  Diagonals<T> cm(numCols);

  // left column
  if (inner_bc[0].first) {
    cm.center[0] = 1;
  } else {
    switch (bcLeft[rowIndex].getType()) {
    case BC_TYPE_CONSTANT: {
      auto [centerCoeffTop, rightCoeffTop] =
          calcBoundaryCoeffConstant(alpha(rowIndex, 0), alpha(rowIndex, 1), sx);
      cm.center[0] = centerCoeffTop;
      cm.right[0] = rightCoeffTop;
      break;
    }
    case BC_TYPE_CLOSED: {
      auto [centerCoeffTop, rightCoeffTop] =
          calcBoundaryCoeffClosed(alpha(rowIndex, 0), alpha(rowIndex, 1), sx);
      cm.center[0] = centerCoeffTop;
      cm.right[0] = rightCoeffTop;
      break;
    }
    default: {
      throw_invalid_argument(
          "Undefined Boundary Condition Type somewhere on Left or Top!");
    }
    }
  }

  // inner columns
  int n = numCols - 1;
  for (int i = 1; i < n; i++) {
    if (inner_bc[i].first) {
      cm.center[i] = 1;
      continue;
    }
    cm.left[i] =
        -sx * calcAlphaIntercell(alpha(rowIndex, i - 1), alpha(rowIndex, i));
    cm.center[i] =
        1 +
        sx * (calcAlphaIntercell(alpha(rowIndex, i), alpha(rowIndex, i + 1)) +
              calcAlphaIntercell(alpha(rowIndex, i - 1), alpha(rowIndex, i)));
    cm.right[i] =
        -sx * calcAlphaIntercell(alpha(rowIndex, i), alpha(rowIndex, i + 1));
  }

  // right column
  if (inner_bc[n].first) {
    cm.center[n] = 1;
  } else {
    switch (bcRight[rowIndex].getType()) {
    case BC_TYPE_CONSTANT: {
      auto [centerCoeffBottom, leftCoeffBottom] = calcBoundaryCoeffConstant(
          alpha(rowIndex, n), alpha(rowIndex, n - 1), sx);
      cm.left[n] = leftCoeffBottom;
      cm.center[n] = centerCoeffBottom;
      break;
    }
    case BC_TYPE_CLOSED: {
      auto [centerCoeffBottom, leftCoeffBottom] = calcBoundaryCoeffClosed(
          alpha(rowIndex, n), alpha(rowIndex, n - 1), sx);
      cm.left[n] = leftCoeffBottom;
      cm.center[n] = centerCoeffBottom;
      break;
    }
    default: {
      throw_invalid_argument(
          "Undefined Boundary Condition Type somewhere on Right or Bottom!");
    }
    }
  }

  return cm;
}

// calculates explicit concentration at boundary in closed case
template <typename T>
constexpr T calcExplicitConcentrationsBoundaryClosed(T conc_center,
                                                     T alpha_center,
                                                     T alpha_neigbor, T sy) {
  return sy * calcAlphaIntercell(alpha_center, alpha_neigbor) * conc_center +
         (1 - sy * (calcAlphaIntercell(alpha_center, alpha_neigbor))) *
             conc_center;
}

// calculates explicity concentration at boundary in constant case
template <typename T>
constexpr T calcExplicitConcentrationsBoundaryConstant(T conc_center, T conc_bc,
                                                       T alpha_center,
                                                       T alpha_neighbor, T sy) {
  const T inter_cell = calcAlphaIntercell(alpha_center, alpha_neighbor);
  return sy * inter_cell * conc_center +
         (1 - sy * (inter_cell + alpha_center)) * conc_center +
         sy * alpha_center * conc_bc;
}

// creates a solution vector for next time step from the current state of
// concentrations
template <class T, class EigenType>
static Eigen::VectorX<T>
createSolutionVector(const EigenType &concentrations,
                     const RowMajMat<T> &alphaX, const RowMajMat<T> &alphaY,
                     const std::vector<BoundaryElement<T>> &bcLeft,
                     const std::vector<BoundaryElement<T>> &bcRight,
                     const std::vector<BoundaryElement<T>> &bcTop,
                     const std::vector<BoundaryElement<T>> &bcBottom,
                     const std::vector<std::pair<bool, T>> &inner_bc,
                     int length, int rowIndex, T sx, T sy) {

  Eigen::VectorX<T> sv(length);
  const std::size_t numRows = concentrations.rows();

  // inner rows
  if (rowIndex > 0 && rowIndex < numRows - 1) {
    for (int i = 0; i < length; i++) {
      if (inner_bc[i].first) {
        sv(i) = inner_bc[i].second;
        continue;
      }
      sv(i) =
          sy *
              calcAlphaIntercell(alphaY(rowIndex, i), alphaY(rowIndex + 1, i)) *
              concentrations(rowIndex + 1, i) +
          (1 - sy * (calcAlphaIntercell(alphaY(rowIndex, i),
                                        alphaY(rowIndex + 1, i)) +
                     calcAlphaIntercell(alphaY(rowIndex - 1, i),
                                        alphaY(rowIndex, i)))) *
              concentrations(rowIndex, i) +
          sy *
              calcAlphaIntercell(alphaY(rowIndex - 1, i), alphaY(rowIndex, i)) *
              concentrations(rowIndex - 1, i);
    }
  }

  // first row
  else if (rowIndex == 0) {
    for (int i = 0; i < length; i++) {
      if (inner_bc[i].first) {
        sv(i) = inner_bc[i].second;
        continue;
      }
      switch (bcTop[i].getType()) {
      case BC_TYPE_CONSTANT: {
        sv(i) = calcExplicitConcentrationsBoundaryConstant(
            concentrations(rowIndex, i), bcTop[i].getValue(),
            alphaY(rowIndex, i), alphaY(rowIndex + 1, i), sy);
        break;
      }
      case BC_TYPE_CLOSED: {
        sv(i) = calcExplicitConcentrationsBoundaryClosed(
            concentrations(rowIndex, i), alphaY(rowIndex, i),
            alphaY(rowIndex + 1, i), sy);
        break;
      }
      default:
        throw_invalid_argument(
            "Undefined Boundary Condition Type somewhere on Left or Top!");
      }
    }
  }

  // last row
  else if (rowIndex == numRows - 1) {
    for (int i = 0; i < length; i++) {
      if (inner_bc[i].first) {
        sv(i) = inner_bc[i].second;
        continue;
      }
      switch (bcBottom[i].getType()) {
      case BC_TYPE_CONSTANT: {
        sv(i) = calcExplicitConcentrationsBoundaryConstant(
            concentrations(rowIndex, i), bcBottom[i].getValue(),
            alphaY(rowIndex, i), alphaY(rowIndex - 1, i), sy);
        break;
      }
      case BC_TYPE_CLOSED: {
        sv(i) = calcExplicitConcentrationsBoundaryClosed(
            concentrations(rowIndex, i), alphaY(rowIndex, i),
            alphaY(rowIndex - 1, i), sy);
        break;
      }
      default:
        throw_invalid_argument(
            "Undefined Boundary Condition Type somewhere on Right or Bottom!");
      }
    }
  }

  // first column -> additional fixed concentration change from perpendicular
  // dimension in constant bc case
  if (bcLeft[rowIndex].getType() == BC_TYPE_CONSTANT && !inner_bc[0].first) {
    sv(0) += 2 * sx * alphaX(rowIndex, 0) * bcLeft[rowIndex].getValue();
  }

  // last column -> additional fixed concentration change from perpendicular
  // dimension in constant bc case
  if (bcRight[rowIndex].getType() == BC_TYPE_CONSTANT &&
      !inner_bc[length - 1].first) {
    sv(length - 1) +=
        2 * sx * alphaX(rowIndex, length - 1) * bcRight[rowIndex].getValue();
  }

  return sv;
}

// solver for linear equation system; A corresponds to coefficient matrix,
// b to the solution vector
// use of EigenLU solver
template <class T>
static Eigen::VectorX<T> EigenLUAlgorithm(Diagonals<T> &A,
                                          Eigen::VectorX<T> &b) {

  // convert A to Eigen sparse matrix
  size_t dim_A = A.center.size();
  Eigen::SparseMatrix<T> A_sparse(dim_A, dim_A);

  A_sparse.insert(0, 0) = A.center[0];
  A_sparse.insert(0, 1) = A.right[0];

  for (size_t i = 1; i < dim_A - 1; i++) {
    A_sparse.insert(i, i - 1) = A.left[i];
    A_sparse.insert(i, i) = A.center[i];
    A_sparse.insert(i, i + 1) = A.right[i];
  }

  A_sparse.insert(dim_A - 1, dim_A - 2) = A.left[dim_A - 1];
  A_sparse.insert(dim_A - 1, dim_A - 1) = A.center[dim_A - 1];

  Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
  solver.analyzePattern(A_sparse);
  solver.factorize(A_sparse);

  return solver.solve(b);
}

// solver for linear equation system; A corresponds to coefficient matrix,
// b to the solution vector
// implementation of Thomas Algorithm
template <class T>
static Eigen::VectorX<T> ThomasAlgorithm(Diagonals<T> &A,
                                         Eigen::VectorX<T> &b) {
  Eigen::Index n = b.size();
  Eigen::VectorX<T> x_vec = b;

  // HACK: write CSV to file
#ifdef WRITE_THOMAS_CSV
#include <fstream>
#include <string>
  static std::uint32_t file_index = 0;
  std::string file_name = "Thomas_" + std::to_string(file_index++) + ".csv";

  std::ofstream out_file;

  out_file.open(file_name, std::ofstream::trunc | std::ofstream::out);

  // print header
  out_file << "Aa, Ab, Ac, b\n";

  // iterate through all elements
  for (std::size_t i = 0; i < n; i++) {
    out_file << a_diag[i] << ", " << b_diag[i] << ", " << c_diag[i] << ", "
             << b[i] << "\n";
  }

  out_file.close();
#endif

  // start solving - c_diag and x_vec are overwritten
  n--;
  A.right[0] /= A.center[0];
  x_vec[0] /= A.center[0];

  for (Eigen::Index i = 1; i < n; i++) {
    A.right[i] /= A.center[i] - A.left[i] * A.right[i - 1];
    x_vec[i] = (x_vec[i] - A.left[i] * x_vec[i - 1]) /
               (A.center[i] - A.left[i] * A.right[i - 1]);
  }

  x_vec[n] = (x_vec[n] - A.left[n] * x_vec[n - 1]) /
             (A.center[n] - A.left[n] * A.right[n - 1]);

  for (Eigen::Index i = n; i-- > 0;) {
    x_vec[i] -= A.right[i] * x_vec[i + 1];
  }

  return x_vec;
}

// BTCS solution for 1D grid
template <class T>
static void BTCS_1D(SimulationInput<T> &input,
                    Eigen::VectorX<T> (*solverFunc)(Diagonals<T> &A,
                                                    Eigen::VectorX<T> &b)) {
  const std::size_t &length = input.colMax;
  T sx = input.timestep / (input.deltaCol * input.deltaCol);

  Eigen::VectorX<T> concentrations_t1(length);

  Diagonals<T> A;
  Eigen::VectorX<T> b(length);

  const auto &alpha = input.alphaX;

  const auto &bc = input.boundaries;

  const auto &bcLeft = bc.getBoundarySide(BC_SIDE_LEFT);
  const auto &bcRight = bc.getBoundarySide(BC_SIDE_RIGHT);

  const auto inner_bc = bc.getInnerBoundaryRow(0);

  RowMajMatMap<T> &concentrations = input.concentrations;
  int rowIndex = 0;
  A = createCoeffMatrix(alpha, bcLeft, bcRight, inner_bc, length, rowIndex,
                        sx); // this is exactly same as in 2D
  for (int i = 0; i < length; i++) {
    b(i) = concentrations(0, i);
  }
  if (bc.getBoundaryElementType(BC_SIDE_LEFT, 0) == BC_TYPE_CONSTANT &&
      !inner_bc[0].first) {
    b(0) += 2 * sx * alpha(0, 0) * bcLeft[0].getValue();
  }
  if (bc.getBoundaryElementType(BC_SIDE_RIGHT, 0) == BC_TYPE_CONSTANT &&
      !inner_bc[length - 1].first) {
    b(length - 1) += 2 * sx * alpha(0, length - 1) * bcRight[0].getValue();
  }

  concentrations = solverFunc(A, b);
}

// BTCS solution for 2D grid
template <class T>
static void BTCS_2D(SimulationInput<T> &input,
                    Eigen::VectorX<T> (*solverFunc)(Diagonals<T> &A,
                                                    Eigen::VectorX<T> &b),
                    int numThreads) {
  const std::size_t &rowMax = input.rowMax;
  const std::size_t &colMax = input.colMax;
  const T sx = input.timestep / (2 * input.deltaCol * input.deltaCol);
  const T sy = input.timestep / (2 * input.deltaRow * input.deltaRow);

  RowMajMat<T> concentrations_t1(rowMax, colMax);

  Diagonals<T> A;
  Eigen::VectorX<T> b;

  const RowMajMat<T> &alphaX = input.alphaX;
  const RowMajMat<T> &alphaY = input.alphaY;

  const auto &bc = input.boundaries;

  const auto &bcLeft = bc.getBoundarySide(BC_SIDE_LEFT);
  const auto &bcRight = bc.getBoundarySide(BC_SIDE_RIGHT);
  const auto &bcTop = bc.getBoundarySide(BC_SIDE_TOP);
  const auto &bcBottom = bc.getBoundarySide(BC_SIDE_BOTTOM);

  RowMajMatMap<T> &concentrations = input.concentrations;

#pragma omp parallel for num_threads(numThreads) private(A, b)
  for (int i = 0; i < rowMax; i++) {
    auto inner_bc = bc.getInnerBoundaryRow(i);

    A = createCoeffMatrix(alphaX, bcLeft, bcRight, inner_bc, colMax, i, sx);
    b = createSolutionVector(concentrations, alphaX, alphaY, bcLeft, bcRight,
                             bcTop, bcBottom, inner_bc, colMax, i, sx, sy);

    concentrations_t1.row(i) = solverFunc(A, b);
  }

  concentrations_t1.transposeInPlace();
  const RowMajMat<T> alphaX_t = alphaX.transpose();
  const RowMajMat<T> alphaY_t = alphaY.transpose();

#pragma omp parallel for num_threads(numThreads) private(A, b)
  for (int i = 0; i < colMax; i++) {
    auto inner_bc = bc.getInnerBoundaryCol(i);
    // swap alphas, boundary conditions and sx/sy for column-wise calculation
    A = createCoeffMatrix(alphaY_t, bcTop, bcBottom, inner_bc, rowMax, i, sy);
    b = createSolutionVector(concentrations_t1, alphaY_t, alphaX_t, bcTop,
                             bcBottom, bcLeft, bcRight, inner_bc, rowMax, i, sy,
                             sx);

    concentrations.col(i) = solverFunc(A, b);
  }
}

// entry point for EigenLU solver; differentiate between 1D and 2D grid
template <class T> void BTCS_LU(SimulationInput<T> &input, int numThreads) {
  tug_assert(input.dim <= 2,
             "Error: Only 1- and 2-dimensional grids are defined!");

  if (input.dim == 1) {
    BTCS_1D(input, EigenLUAlgorithm);
  } else {
    BTCS_2D(input.dim, EigenLUAlgorithm, numThreads);
  }
}

// entry point for Thomas algorithm solver; differentiate 1D and 2D grid
template <class T> void BTCS_Thomas(SimulationInput<T> &input, int numThreads) {
  tug_assert(input.dim <= 2,
             "Error: Only 1- and 2-dimensional grids are defined!");

  if (input.dim == 1) {
    BTCS_1D(input, ThomasAlgorithm);
  } else {
    BTCS_2D(input, ThomasAlgorithm, numThreads);
  }
}
} // namespace tug

#endif // BTCS_H_
