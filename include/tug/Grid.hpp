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
#include <Eigen/Sparse>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <cstddef>

namespace tug {

/**
 * @brief Holds a matrix with concenctration and respective matrix/matrices of
 * alpha coefficients.
 *
 * @tparam T Type to be used for matrices, e.g. double or float
 */
template <class T> class Grid {
public:
  /**
   * @brief Construct a new Grid object.
   *
   * Constructs a new Grid object with given concentrations, defined by an
   * Eigen::Matrix. The domain length is per default the same as the length. The
   * alpha coefficients are set to 1. The dimensions of the grid are determined
   * by the given matrix, which can also be an Eigen::Vector for a 1D-Grid.
   *
   * @param concentrations An Eigen3 MatrixX<T> holding the concentrations.
   */
  Grid(const RowMajMat<T> &concentrations) {
    if (concentrations.rows() == 1) {
      this->dim = 1;
      this->domainCol = static_cast<T>(concentrations.cols());
      this->deltaCol = static_cast<T>(this->domainCol) /
                       static_cast<T>(concentrations.cols()); // -> 1

      this->concentrations = concentrations;
      return;
    }

    if (concentrations.cols() == 1) {
      this->dim = 1;
      this->domainCol = static_cast<T>(concentrations.rows());
      this->deltaCol = static_cast<T>(this->domainCol) /
                       static_cast<T>(concentrations.rows()); // -> 1

      this->concentrations = concentrations.transpose();
      return;
    }

    this->dim = 2;
    this->domainRow = static_cast<T>(concentrations.rows());
    this->domainCol = static_cast<T>(concentrations.cols());
    this->deltaRow = static_cast<T>(this->domainRow) /
                     static_cast<T>(concentrations.rows()); // -> 1
    this->deltaCol = static_cast<T>(this->domainCol) /
                     static_cast<T>(concentrations.cols()); // -> 1

    this->concentrations = concentrations;
    // this->alphaX = RowMajMat<T>::Constant(concentrations.rows(),
    //                                       concentrations.cols(),
    //                                       MAT_INIT_VAL);
    // this->alphaY = RowMajMat<T>::Constant(concentrations.rows(),
    //                                       concentrations.cols(),
    //                                       MAT_INIT_VAL);
  }

  /**
   * @brief Construct a new Grid object.
   *
   * Constructs a new 1D Grid object with given concentrations, defined by a
   * pointer to consecutive memory and the length of the array. The domain
   * length is per default the same as the count of grid cells (length of
   * array). The memory region is mapped internally, changes will affect the
   * original array and the memory shall not be freed. There is no check for
   * correct memory size!
   *
   * @param concentrations Pointer to consecutive memory holding concentrations.
   * @param length Length of the array/the 1D grid.
   */
  Grid(T *concentrations, std::size_t length) : dim(1) {
    this->domainCol = static_cast<T>(length); // -> 1
    this->deltaCol =
        static_cast<T>(this->domainCol) / static_cast<T>(length); // -> 1

    this->concentrations = RowMajMatMap<T>(concentrations, 1, length);
  }

  /**
   * @brief Construct a new Grid object.
   *
   * Constructs a new 2D Grid object with given concentrations, defined by a
   * pointer to consecutive memory and the number of rows and columns. The
   * domain size is per default the same as the number of rows and columns. The
   * memory region is mapped internally, changes will affect the original array
   * and the memory shall not be freed. There is no check for correct memory
   * size!
   *
   * @param concentrations Pointer to consecutive memory holding concentrations.
   * @param row Number of rows.
   * @param col Number of columns.
   */
  Grid(T *concentrations, std::size_t row, std::size_t col) : dim(2) {
    this->domainRow = static_cast<T>(row); // -> 1
    this->domainCol = static_cast<T>(col); // -> 1
    this->deltaCol =
        static_cast<T>(this->domainCol) / static_cast<T>(col); // -> 1
    this->deltaRow =
        static_cast<T>(this->domainRow) / static_cast<T>(row); // -> 1

    this->concentrations = RowMajMatMap<T>(concentrations, row, col);
  }

  /**
   * @brief Gets the concentrations matrix for a Grid.
   *
   * @return An Eigen3 matrix holding the concentrations and having
   * the same dimensions as the grid.
   */
  auto &getConcentrations() { return this->concentrations; }

  void initAlpha() {
    this->alphaX = RowMajMat<T>::Constant(
        this->concentrations.rows(), this->concentrations.cols(), MAT_INIT_VAL);
    if (dim > 1) {

      this->alphaY =
          RowMajMat<T>::Constant(this->concentrations.rows(),
                                 this->concentrations.cols(), MAT_INIT_VAL);
    }
  }

  /**
   * @brief Set the alpha coefficients of a 1D-Grid. Grid must be one
   * dimensional.
   *
   * @param alpha An Eigen3 MatrixX<T> with 1 row holding the alpha
   * coefficients. Matrix columns must have same size as length of grid.
   */
  void setAlpha(const RowMajMat<T> &alpha) {
    tug_assert(dim == 1,
               "Grid is not one dimensional, use 2D setter function!");

    tug_assert(
        alpha.rows() == 1 || alpha.cols() == this->concentrations.cols(),
        "Given matrix of alpha coefficients mismatch with Grid dimensions!");

    this->alphaX = alpha;
  }

  /**
   * @brief Set the alpha coefficients of a 1D-Grid. Grid must be one
   * dimensional.
   *
   * @param alpha A pointer to an array holding the alpha coefficients. Array
   * must have correct dimensions as defined in length. There is no check for
   * correct dimensions, so be careful!
   */
  void setAlpha(T *alpha) {
    tug_assert(dim == 1,
               "Grid is not one dimensional, use 2D setter function!");

    RowMajMatMap<T> map(alpha, 1, this->col);
    this->alphaX = map;
  }

  /**
   * @brief Set the alpha coefficients of a 2D-Grid. Grid must be two
   * dimensional.
   *
   * @param alphaX An Eigen3 MatrixX<T> holding the alpha coefficients in
   * x-direction. Matrix must be of same size as the grid.
   * @param alphaY An Eigen3 MatrixX<T> holding the alpha coefficients in
   * y-direction. Matrix must be of same size as the grid.
   */
  void setAlpha(const RowMajMat<T> &alphaX, const RowMajMat<T> &alphaY) {
    tug_assert(dim == 2,
               "Grid is not two dimensional, use 1D setter function!");

    tug_assert(alphaX.rows() == this->concentrations.rows(),
               "Alpha in x-direction "
               "has wrong number of rows!");
    tug_assert(alphaX.cols() == this->concentrations.cols(),
               "Alpha in x-direction has wrong number of columns!");

    tug_assert(alphaY.rows() == this->concentrations.rows(),
               "Alpha in y-direction "
               "has wrong number of rows!");

    tug_assert(alphaY.cols() == this->concentrations.cols(),
               "Alpha in y-direction has wrong number of columns!");

    this->alphaX = alphaX;
    this->alphaY = alphaY;
  }

  /**
   * @brief Set the alpha coefficients of a 2D-Grid. Grid must be two
   * dimensional.
   *
   * @param alphaX A pointer to an array holding the alpha coefficients in
   * x-direction. Array must have correct dimensions as defined in row and col.
   * There is no check for correct dimensions, so be careful!
   * @param alphaY A pointer to an array holding the alpha coefficients in
   * y-direction. Array must have correct dimensions as defined in row and col.
   * There is no check for correct dimensions, so be careful!
   */
  void setAlpha(T *alphaX, T *alphaY) {
    tug_assert(dim == 2, "Grid is not two dimensional, there is no alphaY!");

    RowMajMatMap<T> mapX(alphaX, this->row, this->col);
    RowMajMatMap<T> mapY(alphaY, this->row, this->col);
    this->alphaX = mapX;
    this->alphaY = mapY;
  }

  /**
   * @brief Gets the matrix of alpha coefficients in x-direction of a 2D-Grid.
   * Grid must be two dimensional.
   *
   * @return A matrix holding the alpha coefficients in x-direction.
   */
  const auto &getAlphaX() const {
    tug_assert(this->alphaX.size() > 0, "AlphaX is empty!");
    return this->alphaX;
  }

  /**
   * @brief Gets the matrix of alpha coefficients in y-direction of a 2D-Grid.
   * Grid must be two dimensional.
   *
   * @return A matrix holding the alpha coefficients in y-direction.
   */
  const auto &getAlphaY() const {
    tug_assert(dim == 2, "Grid is not two dimensional, there is no alphaY!");
    tug_assert(this->alphaY.size() > 0, "AlphaY is empty!");

    return this->alphaY;
  }

  /**
   * @brief Gets the dimensions of the grid.
   *
   * @return Dimensions, either 1 or 2.
   */
  int getDim() const { return this->dim; }

  /**
   * @brief Gets the number of rows of the grid.
   *
   * @return Number of rows.
   */
  int getRow() const { return this->concentrations.rows(); }

  /**
   * @brief Gets the number of columns of the grid.
   *
   * @return Number of columns.
   */
  int getCol() const { return this->concentrations.cols(); }

  /**
   * @brief Sets the domain length of a 1D-Grid. Grid must be one dimensional.
   *
   * @param domainLength A double value of the domain length. Must be positive.
   */
  void setDomain(double domainLength) {
    tug_assert(dim == 1, "Grid is not one dimensional, use 2D domain setter!");
    tug_assert(domainLength > 0, "Given domain length is not positive!");

    this->domainCol = domainLength;
    this->deltaCol =
        double(this->domainCol) / double(this->concentrations.cols());
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
    tug_assert(dim == 2, "Grid is not two dimensional, use 1D domain setter!");
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
    tug_assert(dim == 2, "Grid is not two dimensional, there is no delta in "
                         "y-direction!");

    return this->deltaRow;
  }

private:
  int dim;        // 1D or 2D
  T domainCol;    // number of domain columns
  T domainRow{0}; // number of domain rows
  T deltaCol;     // delta in x-direction (between columns)
  T deltaRow;     // delta in y-direction (between rows)

  RowMajMat<T> concentrations; // Matrix holding grid concentrations
  RowMajMat<T> alphaX; // Matrix holding alpha coefficients in x-direction
  RowMajMat<T> alphaY; // Matrix holding alpha coefficients in y-direction

  static constexpr T MAT_INIT_VAL = 0;
};

using Grid64 = Grid<double>;
using Grid32 = Grid<float>;
} // namespace tug
