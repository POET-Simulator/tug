#pragma once

#include <Eigen/Core>

namespace tug {
/**
 * @brief Alias template for a row-major matrix using Eigen library.
 *
 * This alias template defines a type `RowMajMat` which represents a row-major
 * matrix using the Eigen library. It is a template that takes a type `T` as its
 * template parameter. The matrix is dynamically sized with `Eigen::Dynamic` for
 * both rows and columns. The matrix is stored in row-major order.
 *
 * @tparam T The type of the matrix elements.
 */
template <typename T>
using RowMajMat =
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T> using RowMajMatMap = Eigen::Map<RowMajMat<T>>;
} // namespace tug