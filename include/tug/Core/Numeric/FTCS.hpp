/**
 * @file FTCS.hpp
 * @brief Implementation of heterogenous FTCS (forward time-centered space)
 * solution of diffusion equation in 1D and 2D space.
 *
 */

#ifndef FTCS_H_
#define FTCS_H_

#include "../TugUtils.hpp"
#include "tug/Core/Matrix.hpp"

#include <cstring>
#include <tug/Boundary.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

namespace tug {

// calculates horizontal change on one cell independent of boundary type
template <class T>
static inline T calcHorizontalChange(UnfiormGrid<T> &grid, int &row, int &col) {

  return calcAlphaIntercell(grid.getAlphaX()(row, col + 1),
                            grid.getAlphaX()(row, col)) *
             grid.getConcentrations()(row, col + 1) -
         (calcAlphaIntercell(grid.getAlphaX()(row, col + 1),
                             grid.getAlphaX()(row, col)) +
          calcAlphaIntercell(grid.getAlphaX()(row, col - 1),
                             grid.getAlphaX()(row, col))) *
             grid.getConcentrations()(row, col) +
         calcAlphaIntercell(grid.getAlphaX()(row, col - 1),
                            grid.getAlphaX()(row, col)) *
             grid.getConcentrations()(row, col - 1);
}

// calculates vertical change on one cell independent of boundary type
template <class T>
static inline T calcVerticalChange(UnfiormGrid<T> &grid, int &row, int &col) {

  return calcAlphaIntercell(grid.getAlphaY()(row + 1, col),
                            grid.getAlphaY()(row, col)) *
             grid.getConcentrations()(row + 1, col) -
         (calcAlphaIntercell(grid.getAlphaY()(row + 1, col),
                             grid.getAlphaY()(row, col)) +
          calcAlphaIntercell(grid.getAlphaY()(row - 1, col),
                             grid.getAlphaY()(row, col))) *
             grid.getConcentrations()(row, col) +
         calcAlphaIntercell(grid.getAlphaY()(row - 1, col),
                            grid.getAlphaY()(row, col)) *
             grid.getConcentrations()(row - 1, col);
}

// calculates horizontal change on one cell with a constant left boundary
template <class T>
static inline T calcHorizontalChangeLeftBoundaryConstant(UnfiormGrid<T> &grid,
                                                         Boundary<T> &bc,
                                                         int &row, int &col) {

  return calcAlphaIntercell(grid.getAlphaX()(row, col + 1),
                            grid.getAlphaX()(row, col)) *
             grid.getConcentrations()(row, col + 1) -
         (calcAlphaIntercell(grid.getAlphaX()(row, col + 1),
                             grid.getAlphaX()(row, col)) +
          2 * grid.getAlphaX()(row, col)) *
             grid.getConcentrations()(row, col) +
         2 * grid.getAlphaX()(row, col) *
             bc.getBoundaryElementValue(BC_SIDE_LEFT, row);
}

// calculates horizontal change on one cell with a closed left boundary
template <class T>
static inline T calcHorizontalChangeLeftBoundaryClosed(UnfiormGrid<T> &grid,
                                                       int &row, int &col) {

  return calcAlphaIntercell(grid.getAlphaX()(row, col + 1),
                            grid.getAlphaX()(row, col)) *
         (grid.getConcentrations()(row, col + 1) -
          grid.getConcentrations()(row, col));
}

// checks boundary condition type for a cell on the left edge of grid
template <class T>
static inline T calcHorizontalChangeLeftBoundary(UnfiormGrid<T> &grid,
                                                 Boundary<T> &bc, int &row,
                                                 int &col) {
  if (bc.getBoundaryElementType(BC_SIDE_LEFT, row) == BC_TYPE_CONSTANT) {
    return calcHorizontalChangeLeftBoundaryConstant(grid, bc, row, col);
  } else if (bc.getBoundaryElementType(BC_SIDE_LEFT, row) == BC_TYPE_CLOSED) {
    return calcHorizontalChangeLeftBoundaryClosed(grid, row, col);
  } else {
    throw_invalid_argument("Undefined Boundary Condition Type!");
  }
}

// calculates horizontal change on one cell with a constant right boundary
template <class T>
static inline T calcHorizontalChangeRightBoundaryConstant(UnfiormGrid<T> &grid,
                                                          Boundary<T> &bc,
                                                          int &row, int &col) {

  return 2 * grid.getAlphaX()(row, col) *
             bc.getBoundaryElementValue(BC_SIDE_RIGHT, row) -
         (calcAlphaIntercell(grid.getAlphaX()(row, col - 1),
                             grid.getAlphaX()(row, col)) +
          2 * grid.getAlphaX()(row, col)) *
             grid.getConcentrations()(row, col) +
         calcAlphaIntercell(grid.getAlphaX()(row, col - 1),
                            grid.getAlphaX()(row, col)) *
             grid.getConcentrations()(row, col - 1);
}

// calculates horizontal change on one cell with a closed right boundary
template <class T>
static inline T calcHorizontalChangeRightBoundaryClosed(UnfiormGrid<T> &grid,
                                                        int &row, int &col) {

  return -(calcAlphaIntercell(grid.getAlphaX()(row, col - 1),
                              grid.getAlphaX()(row, col)) *
           (grid.getConcentrations()(row, col) -
            grid.getConcentrations()(row, col - 1)));
}

// checks boundary condition type for a cell on the right edge of grid
template <class T>
static inline T calcHorizontalChangeRightBoundary(UnfiormGrid<T> &grid,
                                                  Boundary<T> &bc, int &row,
                                                  int &col) {
  if (bc.getBoundaryElementType(BC_SIDE_RIGHT, row) == BC_TYPE_CONSTANT) {
    return calcHorizontalChangeRightBoundaryConstant(grid, bc, row, col);
  } else if (bc.getBoundaryElementType(BC_SIDE_RIGHT, row) == BC_TYPE_CLOSED) {
    return calcHorizontalChangeRightBoundaryClosed(grid, row, col);
  } else {
    throw_invalid_argument("Undefined Boundary Condition Type!");
  }
}

// calculates vertical change on one cell with a constant top boundary
template <class T>
static inline T calcVerticalChangeTopBoundaryConstant(UnfiormGrid<T> &grid,
                                                      Boundary<T> &bc, int &row,
                                                      int &col) {

  return calcAlphaIntercell(grid.getAlphaY()(row + 1, col),
                            grid.getAlphaY()(row, col)) *
             grid.getConcentrations()(row + 1, col) -
         (calcAlphaIntercell(grid.getAlphaY()(row + 1, col),
                             grid.getAlphaY()(row, col)) +
          2 * grid.getAlphaY()(row, col)) *
             grid.getConcentrations()(row, col) +
         2 * grid.getAlphaY()(row, col) *
             bc.getBoundaryElementValue(BC_SIDE_TOP, col);
}

// calculates vertical change on one cell with a closed top boundary
template <class T>
static inline T calcVerticalChangeTopBoundaryClosed(UnfiormGrid<T> &grid,
                                                    int &row, int &col) {

  return calcAlphaIntercell(grid.getAlphaY()(row + 1, col),
                            grid.getAlphaY()(row, col)) *
         (grid.getConcentrations()(row + 1, col) -
          grid.getConcentrations()(row, col));
}

// checks boundary condition type for a cell on the top edge of grid
template <class T>
static inline T calcVerticalChangeTopBoundary(UnfiormGrid<T> &grid,
                                              Boundary<T> &bc, int &row,
                                              int &col) {
  if (bc.getBoundaryElementType(BC_SIDE_TOP, col) == BC_TYPE_CONSTANT) {
    return calcVerticalChangeTopBoundaryConstant(grid, bc, row, col);
  } else if (bc.getBoundaryElementType(BC_SIDE_TOP, col) == BC_TYPE_CLOSED) {
    return calcVerticalChangeTopBoundaryClosed(grid, row, col);
  } else {
    throw_invalid_argument("Undefined Boundary Condition Type!");
  }
}

// calculates vertical change on one cell with a constant bottom boundary
template <class T>
static inline T calcVerticalChangeBottomBoundaryConstant(UnfiormGrid<T> &grid,
                                                         Boundary<T> &bc,
                                                         int &row, int &col) {

  return 2 * grid.getAlphaY()(row, col) *
             bc.getBoundaryElementValue(BC_SIDE_BOTTOM, col) -
         (calcAlphaIntercell(grid.getAlphaY()(row, col),
                             grid.getAlphaY()(row - 1, col)) +
          2 * grid.getAlphaY()(row, col)) *
             grid.getConcentrations()(row, col) +
         calcAlphaIntercell(grid.getAlphaY()(row, col),
                            grid.getAlphaY()(row - 1, col)) *
             grid.getConcentrations()(row - 1, col);
}

// calculates vertical change on one cell with a closed bottom boundary
template <class T>
static inline T calcVerticalChangeBottomBoundaryClosed(UnfiormGrid<T> &grid,
                                                       int &row, int &col) {

  return -(calcAlphaIntercell(grid.getAlphaY()(row, col),
                              grid.getAlphaY()(row - 1, col)) *
           (grid.getConcentrations()(row, col) -
            grid.getConcentrations()(row - 1, col)));
}

// checks boundary condition type for a cell on the bottom edge of grid
template <class T>
static inline T calcVerticalChangeBottomBoundary(UnfiormGrid<T> &grid,
                                                 Boundary<T> &bc, int &row,
                                                 int &col) {
  if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, col) == BC_TYPE_CONSTANT) {
    return calcVerticalChangeBottomBoundaryConstant(grid, bc, row, col);
  } else if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, col) == BC_TYPE_CLOSED) {
    return calcVerticalChangeBottomBoundaryClosed(grid, row, col);
  } else {
    throw_invalid_argument("Undefined Boundary Condition Type!");
  }
}

// FTCS solution for 1D grid
template <class T>
static void FTCS_1D(UnfiormGrid<T> &grid, Boundary<T> &bc, T timestep) {
  int colMax = grid.getCol();
  T deltaCol = grid.getDeltaCol();

  RowMajMat<T> &concentrations_grid = grid.getConcentrations();
  // matrix for concentrations at time t+1
  RowMajMat<T> concentrations_t1 = RowMajMat<T>::Constant(1, colMax, 0);

  // only one row in 1D case -> row constant at index 0
  int row = 0;

  // inner cells
  // independent of boundary condition type
  for (int col = 1; col < colMax - 1; col++) {
    concentrations_t1(row, col) = concentrations_grid(row, col) +
                                  timestep / (deltaCol * deltaCol) *
                                      (calcHorizontalChange(grid, row, col));
  }

  // left boundary; hold column constant at index 0
  int col = 0;
  concentrations_t1(row, col) =
      concentrations_grid(row, col) +
      timestep / (deltaCol * deltaCol) *
          (calcHorizontalChangeLeftBoundary(grid, bc, row, col));

  // right boundary; hold column constant at max index
  col = colMax - 1;
  concentrations_t1(row, col) =
      concentrations_grid(row, col) +
      timestep / (deltaCol * deltaCol) *
          (calcHorizontalChangeRightBoundary(grid, bc, row, col));

  // overwrite obsolete concentrations
  std::memcpy(concentrations_grid.data(), concentrations_t1.data(),
              colMax * sizeof(T));
}

// FTCS solution for 2D grid
template <class T>
static void FTCS_2D(UnfiormGrid<T> &grid, Boundary<T> &bc, T timestep,
                    int numThreads) {
  int rowMax = grid.getRow();
  int colMax = grid.getCol();
  T deltaRow = grid.getDeltaRow();
  T deltaCol = grid.getDeltaCol();

  RowMajMat<T> &concentrations_grid = grid.getConcentrations();

  // matrix for concentrations at time t+1
  RowMajMat<T> concentrations_t1 = RowMajMat<T>::Constant(rowMax, colMax, 0);

  // inner cells
  // these are independent of the boundary condition type
  // omp_set_num_threads(10);
#pragma omp parallel for num_threads(numThreads)
  for (int row = 1; row < rowMax - 1; row++) {
    for (int col = 1; col < colMax - 1; col++) {
      concentrations_t1(row, col) = concentrations_grid(row, col) +
                                    timestep / (deltaRow * deltaRow) *
                                        (calcVerticalChange(grid, row, col)) +
                                    timestep / (deltaCol * deltaCol) *
                                        (calcHorizontalChange(grid, row, col));
    }
  }

  // boundary conditions
  // left without corners / looping over rows
  // hold column constant at index 0
  int col = 0;
#pragma omp parallel for num_threads(numThreads)
  for (int row = 1; row < rowMax - 1; row++) {
    concentrations_t1(row, col) =
        concentrations_grid(row, col) +
        timestep / (deltaCol * deltaCol) *
            (calcHorizontalChangeLeftBoundary(grid, bc, row, col)) +
        timestep / (deltaRow * deltaRow) * (calcVerticalChange(grid, row, col));
  }

  // right without corners / looping over rows
  // hold column constant at max index
  col = colMax - 1;
#pragma omp parallel for num_threads(numThreads)
  for (int row = 1; row < rowMax - 1; row++) {
    concentrations_t1(row, col) =
        concentrations_grid(row, col) +
        timestep / (deltaCol * deltaCol) *
            (calcHorizontalChangeRightBoundary(grid, bc, row, col)) +
        timestep / (deltaRow * deltaRow) * (calcVerticalChange(grid, row, col));
  }

  // top without corners / looping over columns
  // hold row constant at index 0
  int row = 0;
#pragma omp parallel for num_threads(numThreads)
  for (int col = 1; col < colMax - 1; col++) {
    concentrations_t1(row, col) =
        concentrations_grid(row, col) +
        timestep / (deltaRow * deltaRow) *
            (calcVerticalChangeTopBoundary(grid, bc, row, col)) +
        timestep / (deltaCol * deltaCol) *
            (calcHorizontalChange(grid, row, col));
  }

  // bottom without corners / looping over columns
  // hold row constant at max index
  row = rowMax - 1;
#pragma omp parallel for num_threads(numThreads)
  for (int col = 1; col < colMax - 1; col++) {
    concentrations_t1(row, col) =
        concentrations_grid(row, col) +
        timestep / (deltaRow * deltaRow) *
            (calcVerticalChangeBottomBoundary(grid, bc, row, col)) +
        timestep / (deltaCol * deltaCol) *
            (calcHorizontalChange(grid, row, col));
  }

  // corner top left
  // hold row and column constant at 0
  row = 0;
  col = 0;
  concentrations_t1(row, col) =
      concentrations_grid(row, col) +
      timestep / (deltaCol * deltaCol) *
          (calcHorizontalChangeLeftBoundary(grid, bc, row, col)) +
      timestep / (deltaRow * deltaRow) *
          (calcVerticalChangeTopBoundary(grid, bc, row, col));

  // corner top right
  // hold row constant at 0 and column constant at max index
  row = 0;
  col = colMax - 1;
  concentrations_t1(row, col) =
      concentrations_grid(row, col) +
      timestep / (deltaCol * deltaCol) *
          (calcHorizontalChangeRightBoundary(grid, bc, row, col)) +
      timestep / (deltaRow * deltaRow) *
          (calcVerticalChangeTopBoundary(grid, bc, row, col));

  // corner bottom left
  // hold row constant at max index and column constant at 0
  row = rowMax - 1;
  col = 0;
  concentrations_t1(row, col) =
      concentrations_grid(row, col) +
      timestep / (deltaCol * deltaCol) *
          (calcHorizontalChangeLeftBoundary(grid, bc, row, col)) +
      timestep / (deltaRow * deltaRow) *
          (calcVerticalChangeBottomBoundary(grid, bc, row, col));

  // corner bottom right
  // hold row and column constant at max index
  row = rowMax - 1;
  col = colMax - 1;
  concentrations_t1(row, col) =
      concentrations_grid(row, col) +
      timestep / (deltaCol * deltaCol) *
          (calcHorizontalChangeRightBoundary(grid, bc, row, col)) +
      timestep / (deltaRow * deltaRow) *
          (calcVerticalChangeBottomBoundary(grid, bc, row, col));

  // overwrite obsolete concentrations
  std::memcpy(concentrations_grid.data(), concentrations_t1.data(),
              rowMax * colMax * sizeof(T));
  // }
}

// entry point; differentiate between 1D and 2D grid
template <class T>
void FTCS(UnfiormGrid<T> &grid, Boundary<T> &bc, T timestep, int &numThreads) {
  if (grid.getDim() == 1) {
    FTCS_1D(grid, bc, timestep);
  } else if (grid.getDim() == 2) {
    FTCS_2D(grid, bc, timestep, numThreads);
  } else {
    throw_invalid_argument(
        "Error: Only 1- and 2-dimensional grids are defined!");
  }
}
} // namespace tug

#endif // FTCS_H_
