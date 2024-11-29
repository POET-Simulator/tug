#pragma once

#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>

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
// template <typename T>
// using RowMajMat =
//     Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <typename T>
class RowMajMat
    : public Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> {
protected:
  std::uint8_t dim;

public:
  RowMajMat(Eigen::Index rows, Eigen::Index cols, T initial_value) : dim(2) {
    *this = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajor>::Constant(rows, cols, initial_value);
  };

  RowMajMat(Eigen::Index n_cells, T initial_value) : dim(1) {
    *this = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic,
                          Eigen::RowMajor>::Constant(1, n_cells, initial_value);
  };

  /**
   * @brief Gets the number of rows of the grid.
   *
   * @return Number of rows.
   */
  Eigen::Index getRow() const { return this->rows(); }

  /**
   * @brief Gets the number of columns of the grid.
   *
   * @return Number of columns.
   */
  Eigen::Index getCol() const { return this->cols(); }

  /**
   * @brief Gets the dimensions of the grid.
   *
   * @return Dimensions, either 1 or 2.
   */
  int getDim() const { return this->dim; }
};

template <typename T> using RowMajMatMap = Eigen::Map<RowMajMat<T>>;
} // namespace tug