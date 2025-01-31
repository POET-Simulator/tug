/**
 * @file FTCS.hpp
 * @brief Implementation of heterogenous FTCS (forward time-centered space)
 * solution of diffusion equation in 1D and 2D space.
 *
 */

#ifndef FTCS_H_
#define FTCS_H_

#include "tug/Core/TugUtils.hpp"
#include <cstddef>
#include <cstring>
#include <tug/Boundary.hpp>
#include <tug/Core/Matrix.hpp>
#include <tug/Core/Numeric/SimulationInput.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

namespace tug {

template <class T>
constexpr T calcChangeInner(T conc_c, T conc_left, T conc_right, T alpha_c,
                            T alpha_left, T alpha_right) {
  const T alpha_center_left = calcAlphaIntercell(alpha_left, alpha_c);
  const T alpha_center_right = calcAlphaIntercell(alpha_right, alpha_c);

  return alpha_center_left * conc_left -
         (alpha_center_left + alpha_center_right) * conc_c +
         alpha_center_right * conc_right;
}

// calculates horizontal change on one cell independent of boundary type
// template <class T>
// static inline T calcHorizontalChange(const RowMajMatMap<T> &concentrations,
//                                      const RowMajMatMap<T> &alphaX, int &row,
//                                      int &col) {
//   return calcAlphaIntercell(alphaX(row, col + 1), alphaX(row, col)) *
//              concentrations(row, col + 1) -
//          (calcAlphaIntercell(alphaX(row, col + 1), alphaX(row, col)) +
//           calcAlphaIntercell(alphaX(row, col - 1), alphaX(row, col))) *
//              concentrations(row, col) +
//          calcAlphaIntercell(alphaX(row, col - 1), alphaX(row, col)) *
//              concentrations(row, col - 1);
// }
//
// // calculates vertical change on one cell independent of boundary type
// template <class T>
// static inline T calcVerticalChange(const RowMajMatMap<T> &concentrations,
//                                    const RowMajMatMap<T> &alphaY, int &row,
//                                    int &col) {
//   return calcAlphaIntercell(alphaY(row + 1, col), alphaY(row, col)) *
//              concentrations(row + 1, col) -
//          (calcAlphaIntercell(alphaY(row + 1, col), alphaY(row, col)) +
//           calcAlphaIntercell(alphaY(row - 1, col), alphaY(row, col))) *
//              concentrations(row, col) +
//          calcAlphaIntercell(alphaY(row - 1, col), alphaY(row, col)) *
//              concentrations(row - 1, col);
// }

// calculates horizontal change on one cell with a constant left boundary
// template <class T>
// static inline T calcHorizontalChangeLeftBoundaryConstant(
//     const RowMajMatMap<T> &concentrations, const RowMajMatMap<T> &alphaX,
//     const Boundary<T> &bc, int &row, int &col) {
//   return calcAlphaIntercell(alphaX(row, col + 1), alphaX(row, col)) *
//              concentrations(row, col + 1) -
//          (calcAlphaIntercell(alphaX(row, col + 1), alphaX(row, col)) +
//           2 * alphaX(row, col)) *
//              concentrations(row, col) +
//          2 * alphaX(row, col) * bc.getBoundaryElementValue(BC_SIDE_LEFT,
//          row);
// }
//
// // calculates horizontal change on one cell with a closed left boundary
// template <class T>
// static inline T
// calcHorizontalChangeLeftBoundaryClosed(const RowMajMatMap<T> &concentrations,
//                                        const RowMajMatMap<T> &alphaX, int
//                                        &row, int &col) {
//   return calcAlphaIntercell(alphaX(row, col + 1), alphaX(row, col)) *
//          (concentrations(row, col + 1) - concentrations(row, col));
// }
//
// // checks boundary condition type for a cell on the left edge of grid
// template <class T>
// static inline T
// calcHorizontalChangeLeftBoundary(const RowMajMatMap<T> &concentrations,
//                                  const RowMajMatMap<T> &alphaX, Boundary<T>
//                                  &bc, int &row, int &col) {
//   if (bc.getBoundaryElementType(BC_SIDE_LEFT, row) == BC_TYPE_CONSTANT) {
//     return calcHorizontalChangeLeftBoundaryConstant(concentrations, alphaX,
//     bc,
//                                                     row, col);
//   } else if (bc.getBoundaryElementType(BC_SIDE_LEFT, row) == BC_TYPE_CLOSED)
//   {
//     return calcHorizontalChangeLeftBoundaryClosed(concentrations, alphaX,
//     row,
//                                                   col);
//   } else {
//     throw_invalid_argument("Undefined Boundary Condition Type!");
//   }
// }
//
// // calculates horizontal change on one cell with a constant right boundary
// template <class T>
// static inline T
// calcHorizontalChangeRightBoundaryConstant(const RowMajMatMap<T>
// &concentrations,
//                                           const RowMajMatMap<T> &alphaX,
//                                           Boundary<T> &bc, int &row, int
//                                           &col) {
//   return 2 * alphaX(row, col) * bc.getBoundaryElementValue(BC_SIDE_RIGHT,
//   row) -
//          (calcAlphaIntercell(alphaX(row, col - 1), alphaX(row, col)) +
//           2 * alphaX(row, col)) *
//              concentrations(row, col) +
//          calcAlphaIntercell(alphaX(row, col - 1), alphaX(row, col)) *
//              concentrations(row, col - 1);
// }
//
// // calculates horizontal change on one cell with a closed right boundary
// template <class T>
// static inline T
// calcHorizontalChangeRightBoundaryClosed(const RowMajMatMap<T>
// &concentrations,
//                                         const RowMajMatMap<T> &alphaX, int
//                                         &row, int &col) {
//   return -(calcAlphaIntercell(alphaX(row, col - 1), alphaX(row, col)) *
//            (concentrations(row, col) - concentrations(row, col - 1)));
// }
//
// // checks boundary condition type for a cell on the right edge of grid
// template <class T>
// static inline T
// calcHorizontalChangeRightBoundary(const RowMajMatMap<T> &concentrations,
//                                   const RowMajMatMap<T> &alphaX,
//                                   Boundary<T> &bc, int &row, int &col) {
//   if (bc.getBoundaryElementType(BC_SIDE_RIGHT, row) == BC_TYPE_CONSTANT) {
//     return calcHorizontalChangeRightBoundaryConstant(concentrations, alphaX,
//     bc,
//                                                      row, col);
//   } else if (bc.getBoundaryElementType(BC_SIDE_RIGHT, row) == BC_TYPE_CLOSED)
//   {
//     return calcHorizontalChangeRightBoundaryClosed(concentrations, alphaX,
//     row,
//                                                    col);
//   } else {
//     throw_invalid_argument("Undefined Boundary Condition Type!");
//   }
// }
//
// // calculates vertical change on one cell with a constant top boundary
// template <class T>
// static inline T
// calcVerticalChangeTopBoundaryConstant(const RowMajMatMap<T> &concentrations,
//                                       const RowMajMatMap<T> &alphaY,
//                                       Boundary<T> &bc, int &row, int &col) {
//   return calcAlphaIntercell(alphaY(row + 1, col), alphaY(row, col)) *
//              concentrations(row + 1, col) -
//          (calcAlphaIntercell(alphaY(row + 1, col), alphaY(row, col)) +
//           2 * alphaY(row, col)) *
//              concentrations(row, col) +
//          2 * alphaY(row, col) * bc.getBoundaryElementValue(BC_SIDE_TOP, col);
// }
//
// // calculates vertical change on one cell with a closed top boundary
// template <class T>
// static inline T
// calcVerticalChangeTopBoundaryClosed(const RowMajMatMap<T> &concentrations,
//                                     const RowMajMatMap<T> &alphaY, int &row,
//                                     int &col) {
//   return calcAlphaIntercell(alphaY(row + 1, col), alphaY(row, col)) *
//          (concentrations(row + 1, col) - concentrations(row, col));
// }
//
// // checks boundary condition type for a cell on the top edge of grid
// template <class T>
// static inline T
// calcVerticalChangeTopBoundary(const RowMajMatMap<T> &concentrations,
//                               const RowMajMatMap<T> &alphaY, Boundary<T> &bc,
//                               int &row, int &col) {
//   if (bc.getBoundaryElementType(BC_SIDE_TOP, col) == BC_TYPE_CONSTANT) {
//     return calcVerticalChangeTopBoundaryConstant(concentrations, alphaY, bc,
//                                                  row, col);
//   } else if (bc.getBoundaryElementType(BC_SIDE_TOP, col) == BC_TYPE_CLOSED) {
//     return calcVerticalChangeTopBoundaryClosed(concentrations, alphaY, row,
//                                                col);
//   } else {
//     throw_invalid_argument("Undefined Boundary Condition Type!");
//   }
// }
//
// // calculates vertical change on one cell with a constant bottom boundary
// template <class T>
// static inline T
// calcVerticalChangeBottomBoundaryConstant(const RowMajMatMap<T>
// &concentrations,
//                                          const RowMajMatMap<T> &alphaY,
//                                          Boundary<T> &bc, int &row, int &col)
//                                          {
//   return 2 * alphaY(row, col) *
//              bc.getBoundaryElementValue(BC_SIDE_BOTTOM, col) -
//          (calcAlphaIntercell(alphaY(row, col), alphaY(row - 1, col)) +
//           2 * alphaY(row, col)) *
//              concentrations(row, col) +
//          calcAlphaIntercell(alphaY(row, col), alphaY(row - 1, col)) *
//              concentrations(row - 1, col);
// }
//
// // calculates vertical change on one cell with a closed bottom boundary
// template <class T>
// static inline T
// calcVerticalChangeBottomBoundaryClosed(const RowMajMatMap<T> &concentrations,
//                                        const RowMajMatMap<T> &alphaY, int
//                                        &row, int &col) {
//   return -(calcAlphaIntercell(alphaY(row, col), alphaY(row - 1, col)) *
//            (concentrations(row, col) - concentrations(row - 1, col)));
// }
//
// // checks boundary condition type for a cell on the bottom edge of grid
// template <class T>
// static inline T
// calcVerticalChangeBottomBoundary(const RowMajMatMap<T> &concentrations,
//                                  const RowMajMatMap<T> &alphaY, Boundary<T>
//                                  &bc, int &row, int &col) {
//   if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, col) == BC_TYPE_CONSTANT) {
//     return calcVerticalChangeBottomBoundaryConstant(concentrations, alphaY,
//     bc,
//                                                     row, col);
//   } else if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, col) ==
//   BC_TYPE_CLOSED) {
//     return calcVerticalChangeBottomBoundaryClosed(concentrations, alphaY,
//     row,
//                                                   col);
//   } else {
//     throw_invalid_argument("Undefined Boundary Condition Type!");
//   }
// }

template <class T>
constexpr T calcChangeBoundary(T conc_c, T conc_neighbor, T alpha_center,
                               T alpha_neighbor, const BoundaryElement<T> &bc) {
  const T alpha_center_neighbor =
      calcAlphaIntercell(alpha_center, alpha_neighbor);
  const T &conc_boundary = bc.getValue();

  switch (bc.getType()) {
  case BC_TYPE_CONSTANT: {
    return 2 * alpha_center * conc_boundary -
           (alpha_center_neighbor + 2 * alpha_center) * conc_c +
           alpha_center_neighbor * conc_neighbor;
  }
  case BC_TYPE_CLOSED: {
    return (alpha_center_neighbor * (conc_neighbor - conc_c));
  }
  }

  tug_assert(false, "Undefined Boundary Condition Type!");
}

// FTCS solution for 1D grid
template <class T> static void FTCS_1D(SimulationInput<T> &input) {
  const std::size_t &colMax = input.colMax;
  const T &deltaCol = input.deltaCol;
  const T &timestep = input.timestep;

  RowMajMatMap<T> &concentrations_grid = input.concentrations;
  // matrix for concentrations at time t+1
  RowMajMat<T> concentrations_t1 = concentrations_grid;

  const auto &alphaX = input.alphaX;
  const auto &bc = input.boundaries;

  // only one row in 1D case -> row constant at index 0
  int row = 0;

  // inner cells
  // independent of boundary condition type
  for (int col = 1; col < colMax - 1; col++) {
    const T &conc_c = concentrations_grid(row, col);
    const T &conc_left = concentrations_grid(row, col - 1);
    const T &conc_right = concentrations_grid(row, col + 1);

    const T &alpha_c = alphaX(row, col);
    const T &alpha_left = alphaX(row, col - 1);
    const T &alpha_right = alphaX(row, col + 1);

    concentrations_t1(row, col) =
        concentrations_grid(row, col) +
        timestep / (deltaCol * deltaCol) *
            calcChangeInner(conc_c, conc_left, conc_right, alpha_c, alpha_left,
                            alpha_right);
  }

  // left boundary; hold column constant at index 0
  {
    int col = 0;
    const T &conc_c = concentrations_grid(row, col);
    const T &conc_right = concentrations_grid(row, col + 1);
    const T &alpha_c = alphaX(row, col);
    const T &alpha_right = alphaX(row, col + 1);
    const BoundaryElement<T> &bc_element =
        input.boundaries.getBoundaryElement(BC_SIDE_LEFT, row);

    concentrations_t1(row, col) =
        concentrations_grid(row, col) +
        timestep / (deltaCol * deltaCol) *
            calcChangeBoundary(conc_c, conc_right, alpha_c, alpha_right,
                               bc_element);
  }

  // right boundary; hold column constant at max index
  {
    int col = colMax - 1;
    const T &conc_c = concentrations_grid(row, col);
    const T &conc_left = concentrations_grid(row, col - 1);
    const T &alpha_c = alphaX(row, col);
    const T &alpha_left = alphaX(row, col - 1);
    const BoundaryElement<T> &bc_element =
        bc.getBoundaryElement(BC_SIDE_RIGHT, row);

    concentrations_t1(row, col) =
        concentrations_grid(row, col) +
        timestep / (deltaCol * deltaCol) *
            calcChangeBoundary(conc_c, conc_left, alpha_c, alpha_left,
                               bc_element);
  }
  // overwrite obsolete concentrations
  concentrations_grid = concentrations_t1;
}

// FTCS solution for 2D grid
template <class T>
static void FTCS_2D(SimulationInput<T> &input, int numThreads) {
  const std::size_t &rowMax = input.rowMax;
  const std::size_t &colMax = input.colMax;
  const T &deltaRow = input.deltaRow;
  const T &deltaCol = input.deltaCol;
  const T &timestep = input.timestep;

  RowMajMatMap<T> &concentrations_grid = input.concentrations;

  // matrix for concentrations at time t+1
  RowMajMat<T> concentrations_t1 = concentrations_grid;

  const auto &alphaX = input.alphaX;
  const auto &alphaY = input.alphaY;
  const auto &bc = input.boundaries;

  const T sx = timestep / (deltaCol * deltaCol);
  const T sy = timestep / (deltaRow * deltaRow);

#pragma omp parallel for num_threads(numThreads)
  for (std::size_t row_i = 0; row_i < rowMax; row_i++) {
    for (std::size_t col_i = 0; col_i < colMax; col_i++) {
      // horizontal change
      T horizontal_change;
      {

        const T &conc_c = concentrations_grid(row_i, col_i);
        const T &alpha_c = alphaX(row_i, col_i);

        if (col_i == 0 || col_i == colMax - 1) {
          // left or right boundary
          const T &conc_neigbor =
              concentrations_grid(row_i, col_i == 0 ? col_i + 1 : col_i - 1);
          const T &alpha_neigbor =
              alphaX(row_i, col_i == 0 ? col_i + 1 : col_i - 1);

          const BoundaryElement<T> &bc_element = bc.getBoundaryElement(
              col_i == 0 ? BC_SIDE_LEFT : BC_SIDE_RIGHT, row_i);

          horizontal_change = calcChangeBoundary(conc_c, conc_neigbor, alpha_c,
                                                 alpha_neigbor, bc_element);
        } else {
          // inner cell
          const T &conc_left = concentrations_grid(row_i, col_i - 1);
          const T &conc_right = concentrations_grid(row_i, col_i + 1);

          const T &alpha_left = alphaX(row_i, col_i - 1);
          const T &alpha_right = alphaX(row_i, col_i + 1);

          horizontal_change = calcChangeInner(conc_c, conc_left, conc_right,
                                              alpha_c, alpha_left, alpha_right);
        }
      }

      // vertical change
      T vertical_change;
      {
        const T &conc_c = concentrations_grid(row_i, col_i);
        const T &alpha_c = alphaY(row_i, col_i);

        if (row_i == 0 || row_i == rowMax - 1) {
          // top or bottom boundary
          const T &conc_neigbor =
              concentrations_grid(row_i == 0 ? row_i + 1 : row_i - 1, col_i);

          const T &alpha_neigbor =
              alphaY(row_i == 0 ? row_i + 1 : row_i - 1, col_i);

          const BoundaryElement<T> &bc_element = bc.getBoundaryElement(
              row_i == 0 ? BC_SIDE_TOP : BC_SIDE_BOTTOM, col_i);

          vertical_change = calcChangeBoundary(conc_c, conc_neigbor, alpha_c,
                                               alpha_neigbor, bc_element);
        } else {
          // inner cell
          const T &conc_bottom = concentrations_grid(row_i - 1, col_i);
          const T &conc_top = concentrations_grid(row_i + 1, col_i);

          const T &alpha_bottom = alphaY(row_i - 1, col_i);
          const T &alpha_top = alphaY(row_i + 1, col_i);

          vertical_change = calcChangeInner(conc_c, conc_bottom, conc_top,
                                            alpha_c, alpha_bottom, alpha_top);
        }
      }

      concentrations_t1(row_i, col_i) = concentrations_grid(row_i, col_i) +
                                        sx * horizontal_change +
                                        sy * vertical_change;
    }
  }

  // inner cells
  // these are independent of the boundary condition type
  // omp_set_num_threads(10);
  // #pragma omp parallel for num_threads(numThreads)
  //   for (int row = 1; row < rowMax - 1; row++) {
  //     for (int col = 1; col < colMax - 1; col++) {
  //       const T &conc_c = concentrations_grid(row, col);
  //       const T &conc_left = concentrations_grid(row, col - 1);
  //       const T &conc_right = concentrations_grid(row, col + 1);
  //       const T &conc_top = concentrations_grid(row + 1, col);
  //       const T &conc_bottom = concentrations_grid(row - 1, col);
  //
  //       const T &alpha_c = alphaX(row, col);
  //       const T &alpha_left = alphaX(row, col - 1);
  //       const T &alpha_right = alphaX(row, col + 1);
  //       const T &alpha_top = alphaY(row + 1, col);
  //       const T &alpha_bottom = alphaY(row - 1, col);
  //
  //       const T horizontal_change = calcChangesInner(
  //           conc_c, conc_left, conc_right, alpha_c, alpha_left, alpha_right);
  //
  //       const T vertical_change = calcChangesInner(
  //           conc_c, conc_bottom, conc_top, alpha_c, alpha_bottom, alpha_top);
  //
  //       concentrations_t1(row, col) =
  //           concentrations_grid(row, col) +
  //           timestep / (deltaRow * deltaRow) * vertical_change +
  //           timestep / (deltaCol * deltaCol) * horizontal_change;
  //     }
  //   }

  // boundary conditions
  // left without corners / looping over rows
  // hold column constant at index 0
  //   int col = 0;
  // #pragma omp parallel for num_threads(numThreads)
  //   for (int row = 1; row < rowMax - 1; row++) {
  //     const T horizontal_change = calcChangeBoundary(
  //         concentrations_grid(row, col), concentrations_grid(row, col + 1),
  //         alphaX(row, col), alphaX(row, col + 1),
  //         bc.getBoundaryElement(BC_SIDE_LEFT, row));
  //
  //     const T vertical_change = calcChangesInner(
  //         concentrations_grid(row, col), concentrations_grid(row - 1, col),
  //         concentrations_grid(row + 1, col), alphaX(row, col),
  //         alphaY(row - 1, col), alphaY(row + 1, col));
  //
  //     concentrations_t1(row, col) =
  //         concentrations_grid(row, col) +
  //         timestep / (deltaCol * deltaCol) * horizontal_change +
  //         timestep / (deltaRow * deltaRow) * vertical_change;
  //   }
  //
  //   // right without corners / looping over rows
  //   // hold column constant at max index
  //   col = colMax - 1;
  // #pragma omp parallel for num_threads(numThreads)
  //   for (int row = 1; row < rowMax - 1; row++) {
  //
  //     concentrations_t1(row, col) =
  //         concentrations_grid(row, col) +
  //         timestep / (deltaCol * deltaCol) *
  //             (calcHorizontalChangeRightBoundary(concentrations_grid, alphaX,
  //             bc,
  //                                                row, col)) +
  //         timestep / (deltaRow * deltaRow) *
  //             (calcVerticalChange(concentrations_grid, alphaY, row, col));
  //   }
  //
  //   // top without corners / looping over columns
  //   // hold row constant at index 0
  //   int row = 0;
  // #pragma omp parallel for num_threads(numThreads)
  //   for (int col = 1; col < colMax - 1; col++) {
  //     concentrations_t1(row, col) =
  //         concentrations_grid(row, col) +
  //         timestep / (deltaRow * deltaRow) *
  //             (calcVerticalChangeTopBoundary(concentrations_grid, alphaY, bc,
  //             row,
  //                                            col)) +
  //         timestep / (deltaCol * deltaCol) *
  //             (calcHorizontalChange(concentrations_grid, alphaX, row, col));
  //   }
  //
  //   // bottom without corners / looping over columns
  //   // hold row constant at max index
  //   row = rowMax - 1;
  // #pragma omp parallel for num_threads(numThreads)
  //   for (int col = 1; col < colMax - 1; col++) {
  //     concentrations_t1(row, col) =
  //         concentrations_grid(row, col) +
  //         timestep / (deltaRow * deltaRow) *
  //             (calcVerticalChangeBottomBoundary(concentrations_grid, alphaY,
  //             bc,
  //                                               row, col)) +
  //         timestep / (deltaCol * deltaCol) *
  //             (calcHorizontalChange(concentrations_grid, alphaX, row, col));
  //   }
  //
  //   // corner top left
  //   // hold row and column constant at 0
  //   row = 0;
  //   col = 0;
  //   concentrations_t1(row, col) =
  //       concentrations_grid(row, col) +
  //       timestep / (deltaCol * deltaCol) *
  //           (calcHorizontalChangeLeftBoundary(concentrations_grid, alphaX,
  //           bc,
  //                                             row, col)) +
  //       timestep / (deltaRow * deltaRow) *
  //           (calcVerticalChangeTopBoundary(concentrations_grid, alphaY, bc,
  //           row,
  //                                          col));
  //
  //   // corner top right
  //   // hold row constant at 0 and column constant at max index
  //   row = 0;
  //   col = colMax - 1;
  //   concentrations_t1(row, col) =
  //       concentrations_grid(row, col) +
  //       timestep / (deltaCol * deltaCol) *
  //           (calcHorizontalChangeRightBoundary(concentrations_grid, alphaX,
  //           bc,
  //                                              row, col)) +
  //       timestep / (deltaRow * deltaRow) *
  //           (calcVerticalChangeTopBoundary(concentrations_grid, alphaY, bc,
  //           row,
  //                                          col));
  //
  //   // corner bottom left
  //   // hold row constant at max index and column constant at 0
  //   row = rowMax - 1;
  //   col = 0;
  //   concentrations_t1(row, col) =
  //       concentrations_grid(row, col) +
  //       timestep / (deltaCol * deltaCol) *
  //           (calcHorizontalChangeLeftBoundary(concentrations_grid, alphaX,
  //           bc,
  //                                             row, col)) +
  //       timestep / (deltaRow * deltaRow) *
  //           (calcVerticalChangeBottomBoundary(concentrations_grid, alphaY,
  //           bc,
  //                                             row, col));
  //
  //   // corner bottom right
  //   // hold row and column constant at max index
  //   row = rowMax - 1;
  //   col = colMax - 1;
  //   concentrations_t1(row, col) =
  //       concentrations_grid(row, col) +
  //       timestep / (deltaCol * deltaCol) *
  //           (calcHorizontalChangeRightBoundary(concentrations_grid, alphaX,
  //           bc,
  //                                              row, col)) +
  //       timestep / (deltaRow * deltaRow) *
  //           (calcVerticalChangeBottomBoundary(concentrations_grid, alphaY,
  //           bc,
  //                                             row, col));

  // overwrite obsolete concentrations
  concentrations_grid = concentrations_t1;
  // }
}

// entry point; differentiate between 1D and 2D grid
template <class T> void FTCS(SimulationInput<T> &input, int &numThreads) {
  tug_assert(input.dim <= 2,
             "Error: Only 1- and 2-dimensional grids are defined!");

  if (input.dim == 1) {
    FTCS_1D(input);
  } else {
    FTCS_2D(input, numThreads);
  }
}
} // namespace tug

#endif // FTCS_H_
