#pragma once

/**
 * @file Grid.hpp
 * @brief API of Grid class, that holds a matrix with concenctrations and a
 *        respective matrix/matrices of alpha coefficients.
 *
 */

#include "tug/UniformGrid.hpp"
#include <Eigen/Core>

namespace tug {

/**
 * @brief Holds a matrix with concenctration and respective matrix/matrices of
 * alpha coefficients.
 *
 * @tparam T Type to be used for matrices, e.g. double or float
 */
template <class T> class UniformAlpha : public RowMajMat<T> {
public:
  UniformAlpha(const UniformGrid<T> &grid, T init_value = 0)
      : RowMajMat<T>::Constant(grid.rows(), grid.cols(), init_value) {}

  UniformAlpha(const UniformGrid<T> &grid, T *alpha)
      : tug::RowMajMatMap<T>(alpha, grid.rows(), grid.cols()) {}
};

} // namespace tug
