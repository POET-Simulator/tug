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

  // overwrite obsolete concentrations
  concentrations_grid = concentrations_t1;
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
