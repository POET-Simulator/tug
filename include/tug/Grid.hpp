#ifndef GRID_H_
#define GRID_H_

/**
 * @file Grid.hpp
 * @brief API of Grid class, that holds a matrix with concenctrations and a
 *        respective matrix/matrices of alpha coefficients.
 *
 */

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <stdexcept>

namespace tug {

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
    this->deltaCol = double(this->domainCol) / double(this->col); // -> 1

    this->concentrations = Eigen::MatrixXd::Constant(1, col, MAT_INIT_VAL);
    this->alphaX = Eigen::MatrixXd::Constant(1, col, MAT_INIT_VAL);
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
    if (row <= 3 || col <= 3) {
      throw std::invalid_argument(
          "Given grid dimensions too small. Must each be greater than 3.");
    }

    this->dim = 2;
    this->deltaRow = double(this->domainRow) / double(this->row); // -> 1
    this->deltaCol = double(this->domainCol) / double(this->col); // -> 1

    this->concentrations = Eigen::MatrixX<T>::Constant(row, col, MAT_INIT_VAL);
    this->alphaX = Eigen::MatrixX<T>::Constant(row, col, MAT_INIT_VAL);
    this->alphaY = Eigen::MatrixX<T>::Constant(row, col, MAT_INIT_VAL);
  }

  /**
   * @brief Sets the concentrations matrix for a 1D or 2D-Grid.
   *
   * @param concentrations An Eigen3 MatrixX<T> holding the concentrations.
   * Matrix must have correct dimensions as defined in row and col. (Or length,
   * in 1D case).
   */
  void setConcentrations(const Eigen::MatrixX<T> &concentrations) {
    if (concentrations.rows() != this->row ||
        concentrations.cols() != this->col) {
      throw std::invalid_argument(
          "Given matrix of concentrations mismatch with Grid dimensions!");
    }

    this->concentrations = concentrations;
  }

  /**
   * @brief Gets the concentrations matrix for a Grid.
   *
   * @return MatrixX<T> An Eigen3 matrix holding the concentrations and having
   * the same dimensions as the grid.
   */
  const Eigen::MatrixX<T> &getConcentrations() { return this->concentrations; }

  /**
   * @brief Set the alpha coefficients of a 1D-Grid. Grid must be one
   * dimensional.
   *
   * @param alpha An Eigen3 MatrixX<T> with 1 row holding the alpha
   * coefficients. Matrix columns must have same size as length of grid.
   */
  void setAlpha(const Eigen::MatrixX<T> &alpha) {
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
   * @brief Set the alpha coefficients of a 2D-Grid. Grid must be two
   * dimensional.
   *
   * @param alphaX An Eigen3 MatrixX<T> holding the alpha coefficients in
   * x-direction. Matrix must be of same size as the grid.
   * @param alphaY An Eigen3 MatrixX<T> holding the alpha coefficients in
   * y-direction. Matrix must be of same size as the grid.
   */
  void setAlpha(const Eigen::MatrixX<T> &alphaX,
                const Eigen::MatrixX<T> &alphaY) {
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
   * @brief Gets the matrix of alpha coefficients of a 1D-Grid. Grid must be one
   * dimensional.
   *
   * @return MatrixX<T> A matrix with 1 row holding the alpha coefficients.
   */
  const Eigen::MatrixX<T> &getAlpha() const {
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
   * @return MatrixX<T> A matrix holding the alpha coefficients in x-direction.
   */
  const Eigen::MatrixX<T> &getAlphaX() const {

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
   * @return MatrixX<T> A matrix holding the alpha coefficients in y-direction.
   */
  const Eigen::MatrixX<T> &getAlphaY() const {

    if (dim != 2) {
      throw std::invalid_argument(
          "Grid is not two dimensional, you should probably use getAlpha()!");
    }

    return this->alphaY;
  }

  /**
   * @brief Gets the dimensions of the grid.
   *
   * @return int Dimensions, either 1 or 2.
   */
  int getDim() const { return this->dim; }

  /**
   * @brief Gets length of 1D grid. Must be one dimensional grid.
   *
   * @return int Length of 1D grid.
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
   * @return int Number of rows.
   */
  int getRow() const { return this->row; }

  /**
   * @brief Gets the number of columns of the grid.
   *
   * @return int Number of columns.
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
   * @return double Delta value.
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
   * @return double Delta value in x-direction.
   */
  T getDeltaCol() const { return this->deltaCol; }

  /**
   * @brief Gets the delta value in y-direction. Must be two dimensional grid.
   *
   * @return double Delta value in y-direction.
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
  int col;                          // number of grid columns
  int row{1};                       // number of grid rows
  int dim;                          // 1D or 2D
  T domainCol;                      // number of domain columns
  T domainRow{0};                   // number of domain rows
  T deltaCol;                       // delta in x-direction (between columns)
  T deltaRow{0};                    // delta in y-direction (between rows)
  Eigen::MatrixX<T> concentrations; // Matrix holding grid concentrations
  Eigen::MatrixX<T> alphaX; // Matrix holding alpha coefficients in x-direction
  Eigen::MatrixX<T> alphaY; // Matrix holding alpha coefficients in y-direction

  static constexpr double MAT_INIT_VAL = 0;
};

using Grid64 = Grid<double>;
using Grid32 = Grid<float>;
} // namespace tug
#endif // GRID_H_
