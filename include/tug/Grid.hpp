#ifndef GRID_H_
#define GRID_H_

/**
 * @file Grid.hpp
 * @brief API of Grid class, that holds a matrix with concenctrations and a
 *        respective matrix/matrices of alpha coefficients.
 *
 */

#include "Core/Matrix.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <stdexcept>

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
   * @brief Constructs a new 1D-Grid object of a given length, which holds a
   * matrix with concentrations and a respective matrix of alpha coefficients.
   *        The domain length is per default the same as the length. The
   * concentrations are all 20 by default and the alpha coefficients are 1.
   *
   * @param length Length of the 1D-Grid. Must be greater than 3.
   */
  Grid(int length) : col(length), domainCol(length) {
    if (length <= 3) {
      throw std::invalid_argument(
          "Given grid length too small. Must be greater than 3.");
    }

    this->dim = 1;
    this->deltaCol =
        static_cast<T>(this->domainCol) / static_cast<T>(this->col); // -> 1

    this->concentrations = RowMajMat<T>::Constant(1, col, MAT_INIT_VAL);
    this->alphaX = RowMajMat<T>::Constant(1, col, MAT_INIT_VAL);
  }

  /**
   * @brief Constructs a new 2D-Grid object of given dimensions, which holds a
   * matrix with concentrations and the respective matrices of alpha coefficient
   * for each direction. The domain in x- and y-direction is per default equal
   * to the col length and row length, respectively. The concentrations are all
   * 20 by default across the entire grid and the alpha coefficients 1 in both
   * directions.
   *
   * @param row Length of the 2D-Grid in y-direction. Must be greater than 3.
   * @param col Length of the 2D-Grid in x-direction. Must be greater than 3.
   */
  Grid(int _row, int _col)
      : row(_row), col(_col), domainRow(_row), domainCol(_col) {
    if (row <= 1 || col <= 1) {
      throw std::invalid_argument(
          "At least one dimension is 1. Use 1D grid for better results.");
    }

    this->dim = 2;
    this->deltaRow =
        static_cast<T>(this->domainRow) / static_cast<T>(this->row); // -> 1
    this->deltaCol =
        static_cast<T>(this->domainCol) / static_cast<T>(this->col); // -> 1

    this->concentrations = RowMajMat<T>::Constant(row, col, MAT_INIT_VAL);
    this->alphaX = RowMajMat<T>::Constant(row, col, MAT_INIT_VAL);
    this->alphaY = RowMajMat<T>::Constant(row, col, MAT_INIT_VAL);
  }

  /**
   * @brief Sets the concentrations matrix for a 1D or 2D-Grid.
   *
   * @param concentrations An Eigen3 MatrixX<T> holding the concentrations.
   * Matrix must have correct dimensions as defined in row and col. (Or length,
   * in 1D case).
   */
  void setConcentrations(const RowMajMat<T> &concentrations) {
    if (concentrations.rows() != this->row ||
        concentrations.cols() != this->col) {
      throw std::invalid_argument(
          "Given matrix of concentrations mismatch with Grid dimensions!");
    }

    this->concentrations = concentrations;
  }

  /**
   * @brief Sets the concentrations matrix for a 1D or 2D-Grid.
   *
   * @param concentrations A pointer to an array holding the concentrations.
   * Array must have correct dimensions as defined in row and col. (Or length,
   * in 1D case). There is no check for correct dimensions, so be careful!
   */
  void setConcentrations(T *concentrations) {
    tug::RowMajMatMap<T> map(concentrations, this->row, this->col);
    this->concentrations = map;
  }

  /**
   * @brief Gets the concentrations matrix for a Grid.
   *
   * @return An Eigen3 matrix holding the concentrations and having
   * the same dimensions as the grid.
   */
  auto &getConcentrations() { return this->concentrations; }

  /**
   * @brief Set the alpha coefficients of a 1D-Grid. Grid must be one
   * dimensional.
   *
   * @param alpha An Eigen3 MatrixX<T> with 1 row holding the alpha
   * coefficients. Matrix columns must have same size as length of grid.
   */
  void setAlpha(const RowMajMat<T> &alpha) {
    if (dim != 1) {
      throw std::invalid_argument(
          "Grid is not one dimensional, you should probably "
          "use 2D setter function!");
    }
    if (alpha.rows() != 1 || alpha.cols() != this->col) {
      throw std::invalid_argument(
          "Given matrix of alpha coefficients mismatch with Grid dimensions!");
    }

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
    if (dim != 1) {
      throw std::invalid_argument(
          "Grid is not one dimensional, you should probably "
          "use 2D setter function!");
    }
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
    if (dim != 2) {
      throw std::invalid_argument(
          "Grid is not two dimensional, you should probably "
          "use 1D setter function!");
    }
    if (alphaX.rows() != this->row || alphaX.cols() != this->col) {
      throw std::invalid_argument(
          "Given matrix of alpha coefficients in x-direction "
          "mismatch with GRid dimensions!");
    }
    if (alphaY.rows() != this->row || alphaY.cols() != this->col) {
      throw std::invalid_argument(
          "Given matrix of alpha coefficients in y-direction "
          "mismatch with GRid dimensions!");
    }

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
    if (dim != 2) {
      throw std::invalid_argument(
          "Grid is not two dimensional, you should probably "
          "use 1D setter function!");
    }
    RowMajMatMap<T> mapX(alphaX, this->row, this->col);
    RowMajMatMap<T> mapY(alphaY, this->row, this->col);
    this->alphaX = mapX;
    this->alphaY = mapY;
  }

  /**
   * @brief Gets the matrix of alpha coefficients of a 1D-Grid. Grid must be one
   * dimensional.
   *
   * @return A matrix with 1 row holding the alpha coefficients.
   */
  const auto &getAlpha() const {
    if (dim != 1) {
      throw std::invalid_argument(
          "Grid is not one dimensional, you should probably "
          "use either getAlphaX() or getAlphaY()!");
    }

    return this->alphaX;
  }

  /**
   * @brief Gets the matrix of alpha coefficients in x-direction of a 2D-Grid.
   * Grid must be two dimensional.
   *
   * @return A matrix holding the alpha coefficients in x-direction.
   */
  const auto &getAlphaX() const {

    if (dim != 2) {
      throw std::invalid_argument(
          "Grid is not two dimensional, you should probably use getAlpha()!");
    }

    return this->alphaX;
  }

  /**
   * @brief Gets the matrix of alpha coefficients in y-direction of a 2D-Grid.
   * Grid must be two dimensional.
   *
   * @return A matrix holding the alpha coefficients in y-direction.
   */
  const auto &getAlphaY() const {

    if (dim != 2) {
      throw std::invalid_argument(
          "Grid is not two dimensional, you should probably use getAlpha()!");
    }

    return this->alphaY;
  }

  /**
   * @brief Gets the dimensions of the grid.
   *
   * @return Dimensions, either 1 or 2.
   */
  int getDim() const { return this->dim; }

  /**
   * @brief Gets length of 1D grid. Must be one dimensional grid.
   *
   * @return Length of 1D grid.
   */
  int getLength() const {
    if (dim != 1) {
      throw std::invalid_argument(
          "Grid is not one dimensional, you should probably "
          "use getRow() or getCol()!");
    }

    return col;
  }

  /**
   * @brief Gets the number of rows of the grid.
   *
   * @return Number of rows.
   */
  int getRow() const { return this->row; }

  /**
   * @brief Gets the number of columns of the grid.
   *
   * @return Number of columns.
   */
  int getCol() const { return this->col; }

  /**
   * @brief Sets the domain length of a 1D-Grid. Grid must be one dimensional.
   *
   * @param domainLength A double value of the domain length. Must be positive.
   */
  void setDomain(double domainLength) {
    if (dim != 1) {
      throw std::invalid_argument(
          "Grid is not one dimensional, you should probaly "
          "use the 2D domain setter!");
    }
    if (domainLength <= 0) {
      throw std::invalid_argument("Given domain length is not positive!");
    }

    this->domainCol = domainLength;
    this->deltaCol = double(this->domainCol) / double(this->col);
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
    if (dim != 2) {
      throw std::invalid_argument(
          "Grid is not two dimensional, you should probably "
          "use the 1D domain setter!");
    }
    if (domainRow <= 0 || domainCol <= 0) {
      throw std::invalid_argument("Given domain size is not positive!");
    }

    this->domainRow = domainRow;
    this->domainCol = domainCol;
    this->deltaRow = double(this->domainRow) / double(this->row);
    this->deltaCol = double(this->domainCol) / double(this->col);
  }

  /**
   * @brief Gets the delta value for 1D-Grid. Grid must be one dimensional.
   *
   * @return Delta value.
   */
  T getDelta() const {

    if (dim != 1) {
      throw std::invalid_argument(
          "Grid is not one dimensional, you should probably "
          "use the 2D delta getters");
    }

    return this->deltaCol;
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
    if (dim != 2) {
      throw std::invalid_argument(
          "Grid is not two dimensional, meaning there is no "
          "delta in y-direction!");
    }

    return this->deltaRow;
  }

private:
  int col;        // number of grid columns
  int row{1};     // number of grid rows
  int dim;        // 1D or 2D
  T domainCol;    // number of domain columns
  T domainRow{0}; // number of domain rows
  T deltaCol;     // delta in x-direction (between columns)
  T deltaRow{0};  // delta in y-direction (between rows)

  RowMajMat<T> concentrations; // Matrix holding grid concentrations
  RowMajMat<T> alphaX; // Matrix holding alpha coefficients in x-direction
  RowMajMat<T> alphaY; // Matrix holding alpha coefficients in y-direction

  static constexpr T MAT_INIT_VAL = 0;
};

using Grid64 = Grid<double>;
using Grid32 = Grid<float>;
} // namespace tug
#endif // GRID_H_
