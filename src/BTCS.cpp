/**
 * @file BTCSv2.cpp
 * @brief Implementation of heterogenous BTCS (backward time-centered space)
 * solution of diffusion equation in 1D and 2D space. Internally the
 * alternating-direction implicit (ADI) method is used. Version 2, because
 * Version 1 was an implementation for the homogeneous BTCS solution.
 *
 */

#include "Schemes.hpp"
#include "TugUtils.hpp"

#include <Eigen/src/Core/util/Meta.h>
#include <cstddef>
#include <tug/Boundary.hpp>
#include <tug/Grid.hpp>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

namespace tug {

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
static Eigen::SparseMatrix<T>
createCoeffMatrix(const Eigen::MatrixX<T> &alpha,
                  const std::vector<BoundaryElement<T>> &bcLeft,
                  const std::vector<BoundaryElement<T>> &bcRight, int numCols,
                  int rowIndex, T sx) {

  // square matrix of column^2 dimension for the coefficients
  Eigen::SparseMatrix<T> cm(numCols, numCols);
  cm.reserve(Eigen::VectorXi::Constant(numCols, 3));

  // left column

  switch (bcLeft[rowIndex].getType()) {
  case BC_TYPE_CONSTANT: {
    auto [centerCoeffTop, rightCoeffTop] =
        calcBoundaryCoeffConstant(alpha(rowIndex, 0), alpha(rowIndex, 1), sx);
    cm.insert(0, 0) = centerCoeffTop;
    cm.insert(0, 1) = rightCoeffTop;
    break;
  }
  case BC_TYPE_CLOSED: {
    auto [centerCoeffTop, rightCoeffTop] =
        calcBoundaryCoeffClosed(alpha(rowIndex, 0), alpha(rowIndex, 1), sx);
    cm.insert(0, 0) = centerCoeffTop;
    cm.insert(0, 1) = rightCoeffTop;
    break;
  }
  default: {
    throw_invalid_argument(
        "Undefined Boundary Condition Type somewhere on Left or Top!");
  }
  }

  // inner columns
  int n = numCols - 1;
  for (int i = 1; i < n; i++) {
    cm.insert(i, i - 1) =
        -sx * calcAlphaIntercell(alpha(rowIndex, i - 1), alpha(rowIndex, i));
    cm.insert(i, i) =
        1 +
        sx * (calcAlphaIntercell(alpha(rowIndex, i), alpha(rowIndex, i + 1)) +
              calcAlphaIntercell(alpha(rowIndex, i - 1), alpha(rowIndex, i)));
    cm.insert(i, i + 1) =
        -sx * calcAlphaIntercell(alpha(rowIndex, i), alpha(rowIndex, i + 1));
  }

  // right column

  switch (bcRight[rowIndex].getType()) {
  case BC_TYPE_CONSTANT: {
    auto [centerCoeffBottom, leftCoeffBottom] = calcBoundaryCoeffConstant(
        alpha(rowIndex, n), alpha(rowIndex, n - 1), sx);
    cm.insert(n, n - 1) = leftCoeffBottom;
    cm.insert(n, n) = centerCoeffBottom;
    break;
  }
  case BC_TYPE_CLOSED: {
    auto [centerCoeffBottom, leftCoeffBottom] =
        calcBoundaryCoeffClosed(alpha(rowIndex, n), alpha(rowIndex, n - 1), sx);
    cm.insert(n, n - 1) = leftCoeffBottom;
    cm.insert(n, n) = centerCoeffBottom;
    break;
  }
  default: {
    throw_invalid_argument(
        "Undefined Boundary Condition Type somewhere on Right or Bottom!");
  }
  }

  cm.makeCompressed(); // important for Eigen solver

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
  return sy * calcAlphaIntercell(alpha_center, alpha_neighbor) * conc_center +
         (1 - sy * (calcAlphaIntercell(alpha_center, alpha_center) +
                    2 * alpha_center)) *
             conc_center +
         sy * alpha_center * conc_bc;
}

// creates a solution vector for next time step from the current state of
// concentrations
template <class T>
static Eigen::VectorX<T>
createSolutionVector(const Eigen::MatrixX<T> &concentrations,
                     const Eigen::MatrixX<T> &alphaX,
                     const Eigen::MatrixX<T> &alphaY,
                     const std::vector<BoundaryElement<T>> &bcLeft,
                     const std::vector<BoundaryElement<T>> &bcRight,
                     const std::vector<BoundaryElement<T>> &bcTop,
                     const std::vector<BoundaryElement<T>> &bcBottom,
                     int length, int rowIndex, T sx, T sy) {

  Eigen::VectorX<T> sv(length);
  const std::size_t numRows = concentrations.rows();

  // inner rows
  if (rowIndex > 0 && rowIndex < numRows - 1) {
    for (int i = 0; i < length; i++) {
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
  if (bcLeft[rowIndex].getType() == BC_TYPE_CONSTANT) {
    sv(0) += 2 * sx * alphaX(rowIndex, 0) * bcLeft[rowIndex].getValue();
  }

  // last column -> additional fixed concentration change from perpendicular
  // dimension in constant bc case
  if (bcRight[rowIndex].getType() == BC_TYPE_CONSTANT) {
    sv(length - 1) +=
        2 * sx * alphaX(rowIndex, length - 1) * bcRight[rowIndex].getValue();
  }

  return sv;
}

// solver for linear equation system; A corresponds to coefficient matrix,
// b to the solution vector
// use of EigenLU solver
template <class T>
static Eigen::VectorX<T> EigenLUAlgorithm(Eigen::SparseMatrix<T> &A,
                                          Eigen::VectorX<T> &b) {

  Eigen::SparseLU<Eigen::SparseMatrix<T>> solver;
  solver.analyzePattern(A);
  solver.factorize(A);

  return solver.solve(b);
}

// solver for linear equation system; A corresponds to coefficient matrix,
// b to the solution vector
// implementation of Thomas Algorithm
template <class T>
static Eigen::VectorX<T> ThomasAlgorithm(Eigen::SparseMatrix<T> &A,
                                         Eigen::VectorX<T> &b) {
  Eigen::Index n = b.size();

  Eigen::VectorX<T> a_diag(n);
  Eigen::VectorX<T> b_diag(n);
  Eigen::VectorX<T> c_diag(n);
  Eigen::VectorX<T> x_vec = b;

  // Fill diagonals vectors
  b_diag[0] = A.coeff(0, 0);
  c_diag[0] = A.coeff(0, 1);

  for (Eigen::Index i = 1; i < n - 1; i++) {
    a_diag[i] = A.coeff(i, i - 1);
    b_diag[i] = A.coeff(i, i);
    c_diag[i] = A.coeff(i, i + 1);
  }
  a_diag[n - 1] = A.coeff(n - 1, n - 2);
  b_diag[n - 1] = A.coeff(n - 1, n - 1);

  // start solving - c_diag and x_vec are overwritten
  n--;
  c_diag[0] /= b_diag[0];
  x_vec[0] /= b_diag[0];

  for (Eigen::Index i = 1; i < n; i++) {
    c_diag[i] /= b_diag[i] - a_diag[i] * c_diag[i - 1];
    x_vec[i] = (x_vec[i] - a_diag[i] * x_vec[i - 1]) /
               (b_diag[i] - a_diag[i] * c_diag[i - 1]);
  }

  x_vec[n] = (x_vec[n] - a_diag[n] * x_vec[n - 1]) /
             (b_diag[n] - a_diag[n] * c_diag[n - 1]);

  for (Eigen::Index i = n; i-- > 0;) {
    x_vec[i] -= c_diag[i] * x_vec[i + 1];
  }

  return x_vec;
}

// BTCS solution for 1D grid
template <class T>
static void BTCS_1D(Grid<T> &grid, Boundary<T> &bc, T timestep,
                    Eigen::VectorX<T> (*solverFunc)(Eigen::SparseMatrix<T> &A,
                                                    Eigen::VectorX<T> &b)) {
  int length = grid.getLength();
  T sx = timestep / (grid.getDelta() * grid.getDelta());

  Eigen::VectorX<T> concentrations_t1(length);

  Eigen::SparseMatrix<T> A;
  Eigen::VectorX<T> b(length);

  const auto &alpha = grid.getAlpha();

  const auto &bcLeft = bc.getBoundarySide(BC_SIDE_LEFT);
  const auto &bcRight = bc.getBoundarySide(BC_SIDE_RIGHT);

  Eigen::MatrixX<T> concentrations = grid.getConcentrations();
  int rowIndex = 0;
  A = createCoeffMatrix(alpha, bcLeft, bcRight, length, rowIndex,
                        sx); // this is exactly same as in 2D
  for (int i = 0; i < length; i++) {
    b(i) = concentrations(0, i);
  }
  if (bc.getBoundaryElementType(BC_SIDE_LEFT, 0) == BC_TYPE_CONSTANT) {
    b(0) += 2 * sx * alpha(0, 0) * bcLeft[0].getValue();
  }
  if (bc.getBoundaryElementType(BC_SIDE_RIGHT, 0) == BC_TYPE_CONSTANT) {
    b(length - 1) += 2 * sx * alpha(0, length - 1) * bcRight[0].getValue();
  }

  concentrations_t1 = solverFunc(A, b);

  for (int j = 0; j < length; j++) {
    concentrations(0, j) = concentrations_t1(j);
  }

  grid.setConcentrations(concentrations);
}

// BTCS solution for 2D grid
template <class T>
static void BTCS_2D(Grid<T> &grid, Boundary<T> &bc, T timestep,
                    Eigen::VectorX<T> (*solverFunc)(Eigen::SparseMatrix<T> &A,
                                                    Eigen::VectorX<T> &b),
                    int numThreads) {
  int rowMax = grid.getRow();
  int colMax = grid.getCol();
  T sx = timestep / (2 * grid.getDeltaCol() * grid.getDeltaCol());
  T sy = timestep / (2 * grid.getDeltaRow() * grid.getDeltaRow());

  Eigen::MatrixX<T> concentrations_t1 =
      Eigen::MatrixX<T>::Constant(rowMax, colMax, 0);
  Eigen::VectorX<T> row_t1(colMax);

  Eigen::SparseMatrix<T> A;
  Eigen::VectorX<T> b;

  auto alphaX = grid.getAlphaX();
  auto alphaY = grid.getAlphaY();

  const auto &bcLeft = bc.getBoundarySide(BC_SIDE_LEFT);
  const auto &bcRight = bc.getBoundarySide(BC_SIDE_RIGHT);
  const auto &bcTop = bc.getBoundarySide(BC_SIDE_TOP);
  const auto &bcBottom = bc.getBoundarySide(BC_SIDE_BOTTOM);

  Eigen::MatrixX<T> concentrations = grid.getConcentrations();

#pragma omp parallel for num_threads(numThreads) private(A, b, row_t1)
  for (int i = 0; i < rowMax; i++) {

    A = createCoeffMatrix(alphaX, bcLeft, bcRight, colMax, i, sx);
    b = createSolutionVector(concentrations, alphaX, alphaY, bcLeft, bcRight,
                             bcTop, bcBottom, colMax, i, sx, sy);

    row_t1 = solverFunc(A, b);

    concentrations_t1.row(i) = row_t1;
  }

  concentrations_t1.transposeInPlace();
  concentrations.transposeInPlace();
  alphaX.transposeInPlace();
  alphaY.transposeInPlace();

#pragma omp parallel for num_threads(numThreads) private(A, b, row_t1)
  for (int i = 0; i < colMax; i++) {
    // swap alphas, boundary conditions and sx/sy for column-wise calculation
    A = createCoeffMatrix(alphaY, bcTop, bcBottom, rowMax, i, sy);
    b = createSolutionVector(concentrations_t1, alphaY, alphaX, bcTop, bcBottom,
                             bcLeft, bcRight, rowMax, i, sy, sx);

    row_t1 = solverFunc(A, b);

    concentrations.row(i) = row_t1;
  }

  concentrations.transposeInPlace();

  grid.setConcentrations(concentrations);
}

// entry point for EigenLU solver; differentiate between 1D and 2D grid
template <class T>
void BTCS_LU(Grid<T> &grid, Boundary<T> &bc, T timestep, int numThreads) {
  if (grid.getDim() == 1) {
    BTCS_1D(grid, bc, timestep, EigenLUAlgorithm);
  } else if (grid.getDim() == 2) {
    BTCS_2D(grid, bc, timestep, EigenLUAlgorithm, numThreads);
  } else {
    throw_invalid_argument(
        "Error: Only 1- and 2-dimensional grids are defined!");
  }
}

// entry point for Thomas algorithm solver; differentiate 1D and 2D grid
template <class T>
void BTCS_Thomas(Grid<T> &grid, Boundary<T> &bc, T timestep, int numThreads) {
  if (grid.getDim() == 1) {
    BTCS_1D(grid, bc, timestep, ThomasAlgorithm);
  } else if (grid.getDim() == 2) {
    BTCS_2D(grid, bc, timestep, ThomasAlgorithm, numThreads);
  } else {
    throw_invalid_argument(
        "Error: Only 1- and 2-dimensional grids are defined!");
  }
}

template void BTCS_Thomas(Grid<double> &grid, Boundary<double> &bc,
                          double timestep, int numThreads);
template void BTCS_Thomas(Grid<float> &grid, Boundary<float> &bc,
                          float timestep, int numThreads);

template void BTCS_LU(Grid<double> &grid, Boundary<double> &bc, double timestep,
                      int numThreads);
template void BTCS_LU(Grid<float> &grid, Boundary<float> &bc, float timestep,
                      int numThreads);
} // namespace tug
