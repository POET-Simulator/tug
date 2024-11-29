#pragma once

/**
 * @file Grid.hpp
 * @brief API of Grid class, that holds a matrix with concenctrations and a
 *        respective matrix/matrices of alpha coefficients.
 *
 */

#include "Core/Matrix.hpp"
#include "tug/Core/TugUtils.hpp"
#include <Eigen/Core>

namespace tug {

/**
 * @brief Holds a matrix with concenctration and respective matrix/matrices of
 * alpha coefficients.
 *
 * @tparam T Type to be used for matrices, e.g. double or float
 */
template <class T> class UniformGrid : public RowMajMat<T> {
public:
  UniformGrid(Eigen::Index n_cells, T spatial_size, T initial_value)
      : RowMajMat<T>(n_cells, initial_value), domainCol(spatial_size) {
    this->deltaCol = static_cast<T>(this->domainCol) / static_cast<T>(n_cells);
  }

  UniformGrid(Eigen::Index row, Eigen::Index col, T spatial_row_size,
              T spatial_col_size, T initial_value)
      : RowMajMat<T>(row, col, initial_value), domainRow(spatial_row_size),
        domainCol(spatial_col_size) {
    this->deltaRow = static_cast<T>(this->domainRow) / static_cast<T>(row);
    this->deltaCol = static_cast<T>(this->domainCol) / static_cast<T>(col);
  }

  UniformGrid(T *concentrations, Eigen::Index n_cells, T spatial_size)
      : RowMajMatMap<T>(concentrations, 1, n_cells), domainCol(spatial_size) {
    this->deltaCol = static_cast<T>(this->domainCol) / static_cast<T>(n_cells);
  }

  UniformGrid(T *concentrations, Eigen::Index row, Eigen::Index col,
              T spatial_row_size, T spatial_col_size)
      : RowMajMatMap<T>(concentrations, row, col), domainRow(spatial_row_size),
        domainCol(spatial_col_size) {
    this->deltaRow = static_cast<T>(this->domainRow) / static_cast<T>(row);
    this->deltaCol = static_cast<T>(this->domainCol) / static_cast<T>(col);
  }

  /**
   * @brief Sets the domain length of a 1D-Grid. Grid must be one dimensional.
   *
   * @param domainLength A double value of the domain length. Must be positive.
   */
  void setDomain(double domainLength) {
    tug_assert(this->dim == 1,
               "Grid is not one dimensional, use 2D domain setter!");
    tug_assert(domainLength > 0, "Given domain length is not positive!");

    this->domainCol = domainLength;
    this->deltaCol = double(this->domainCol) / double(this->cols());
  }

  /**
   * @brief Sets the domain size of a 2D-Grid. Grid must be two dimensional.
   *
   * @param domainRow A double value of the domain size in y-direction. Must
   * be positive.
   * @param domainCol A double value of the domain size in x-direction. Must
   * be positive.
   */
  void setDomain(double domainRow, double domainCol) {
    tug_assert(this->dim == 2,
               "Grid is not two dimensional, use 1D domain setter!");
    tug_assert(domainCol > 0,
               "Given domain size in x-direction is not positive!");
    tug_assert(domainRow > 0,
               "Given domain size in y-direction is not positive!");

    this->domainRow = domainRow;
    this->domainCol = domainCol;
    this->deltaRow =
        double(this->domainRow) / double(this->concentrations.rows());
    this->deltaCol =
        double(this->domainCol) / double(this->concentrations.cols());
  }

  /**
   * @brief Gets the delta value in x-direction.
   *
   * @return Delta value in x-direction.
   */
  T getDeltaCol() const { return this->deltaCol; }

  /**
   * @brief Gets the delta value in y-direction. Must be two dimensional grid.
   *
   * @return Delta value in y-direction.
   */
  T getDeltaRow() const {
    tug_assert(this->dim == 2,
               "Grid is not two dimensional, there is no delta in "
               "y-direction!");

    return this->deltaRow;
  }

private:
  T domainCol;    // number of domain columns
  T domainRow{0}; // number of domain rows
  T deltaCol;     // delta in x-direction (between columns)
  T deltaRow;     // delta in y-direction (between rows)
};

using UniformGrid64 = UniformGrid<double>;
using UniformGrid32 = UniformGrid<float>;
} // namespace tug
