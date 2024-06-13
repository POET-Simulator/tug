/**
 * @file Boundary.hpp
 * @brief API of Boundary class, that holds all information for each boundary
 * condition at the edges of the diffusion grid.
 *
 */
#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "Grid.hpp"
#include "tug/Core/TugUtils.hpp"

#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace tug {

/**
 * @brief Enum defining the two implemented boundary conditions.
 *
 */
enum BC_TYPE { BC_TYPE_CLOSED, BC_TYPE_CONSTANT };

/**
 * @brief Enum defining all 4 possible sides to a 1D and 2D grid.
 *
 */
enum BC_SIDE { BC_SIDE_LEFT, BC_SIDE_RIGHT, BC_SIDE_TOP, BC_SIDE_BOTTOM };

/**
 * This class defines the boundary conditions of individual boundary elements.
 * These can be flexibly used and combined later in other classes.
 * The class serves as an auxiliary class for structuring the Boundary class.
 *
 * @tparam T Data type of the boundary condition element
 */
template <class T> class BoundaryElement {
public:
  /**
   * @brief Construct a new Boundary Element object for the closed case.
   *        The boundary type is here automatically set to the type
   *        BC_TYPE_CLOSED, where the value takes -1 and does not hold any
   *        physical meaning.
   */
  BoundaryElement(){};

  /**
   * @brief Construct a new Boundary Element object for the constant case.
   *        The boundary type is automatically set to the type
   * BC_TYPE_CONSTANT.
   *
   * @param value Value of the constant concentration to be assumed at the
   *              corresponding boundary element.
   */
  BoundaryElement(T _value) : value(_value), type(BC_TYPE_CONSTANT) {}

  /**
   * @brief Allows changing the boundary type of a corresponding
   *        BoundaryElement object.
   *
   * @param type Type of boundary condition. Either BC_TYPE_CONSTANT or
                 BC_TYPE_CLOSED.
   */
  void setType(BC_TYPE type) { this->type = type; };

  /**
   * @brief Sets the value of a boundary condition for the constant case.
   *
   * @param value Concentration to be considered constant for the
   *              corresponding boundary element.
   */
  void setValue(double value) {
    if (type == BC_TYPE_CLOSED) {
      throw std::invalid_argument(
          "No constant boundary concentrations can be set for closed "
          "boundaries. Please change type first.");
    }
    this->value = value;
  }

  /**
   * @brief Return the type of the boundary condition, i.e. whether the
   *        boundary is considered closed or constant.
   *
   * @return Type of boundary condition, either BC_TYPE_CLOSED or
             BC_TYPE_CONSTANT.
   */
  BC_TYPE getType() const { return this->type; }

  /**
   * @brief Return the concentration value for the constant boundary condition.
   *
   * @return Value of the concentration.
   */
  T getValue() const { return this->value; }

private:
  BC_TYPE type{BC_TYPE_CLOSED};
  T value{-1};
};

/**
 * This class implements the functionality and management of the boundary
 * conditions in the grid to be simulated.
 *
 * @tparam Data type of the boundary condition value
 */
template <class T> class Boundary {
public:
  /**
   * @brief Creates a boundary object for a 1D grid
   *
   * @param length Length of the grid
   */
  Boundary(std::uint32_t length) : Boundary(Grid<T>(length)){};

  /**
   * @brief Creates a boundary object for a 2D grid
   *
   * @param rows Number of rows of the grid
   * @param cols Number of columns of the grid
   */
  Boundary(std::uint32_t rows, std::uint32_t cols)
      : Boundary(Grid<T>(rows, cols)){};

  /**
   * @brief Creates a boundary object based on the passed grid object and
   *        initializes the boundaries as closed.
   *
   * @param grid Grid object on the basis of which the simulation takes place
   *             and from which the dimensions (in 2D case) are taken.
   */
  Boundary(const Grid<T> &grid)
      : dim(grid.getDim()), cols(grid.getCol()), rows(grid.getRow()) {
    if (this->dim == 1) {
      this->boundaries = std::vector<std::vector<BoundaryElement<T>>>(
          2); // in 1D only left and right boundary

      this->boundaries[BC_SIDE_LEFT].push_back(BoundaryElement<T>());
      this->boundaries[BC_SIDE_RIGHT].push_back(BoundaryElement<T>());
    } else if (this->dim == 2) {
      this->boundaries = std::vector<std::vector<BoundaryElement<T>>>(4);

      this->boundaries[BC_SIDE_LEFT] =
          std::vector<BoundaryElement<T>>(this->rows, BoundaryElement<T>());
      this->boundaries[BC_SIDE_RIGHT] =
          std::vector<BoundaryElement<T>>(this->rows, BoundaryElement<T>());
      this->boundaries[BC_SIDE_TOP] =
          std::vector<BoundaryElement<T>>(this->cols, BoundaryElement<T>());
      this->boundaries[BC_SIDE_BOTTOM] =
          std::vector<BoundaryElement<T>>(this->cols, BoundaryElement<T>());
    }
  }

  /**
   * @brief Sets all elements of the specified boundary side to the boundary
   *        condition closed.
   *
   * @param side Side to be set to closed, e.g. BC_SIDE_LEFT.
   */
  void setBoundarySideClosed(BC_SIDE side) {
    tug_assert((this->dim > 1) ||
                   ((side == BC_SIDE_LEFT) || (side == BC_SIDE_RIGHT)),
               "For the "
               "one-dimensional case, only the BC_SIDE_LEFT and BC_SIDE_RIGHT "
               "borders exist.");

    const bool is_vertical = side == BC_SIDE_LEFT || side == BC_SIDE_RIGHT;
    const int n = is_vertical ? this->rows : this->cols;

    this->boundaries[side] =
        std::vector<BoundaryElement<T>>(n, BoundaryElement<T>());
  }

  /**
   * @brief Sets all elements of the specified boundary side to the boundary
   *        condition constant. Thereby the concentration values of the
   *        boundaries are set to the passed value.
   *
   * @param side Side to be set to constant, e.g. BC_SIDE_LEFT.
   * @param value Concentration to be set for all elements of the specified
   * page.
   */
  void setBoundarySideConstant(BC_SIDE side, double value) {
    tug_assert((this->dim > 1) ||
                   ((side == BC_SIDE_LEFT) || (side == BC_SIDE_RIGHT)),
               "For the "
               "one-dimensional case, only the BC_SIDE_LEFT and BC_SIDE_RIGHT "
               "borders exist.");

    const bool is_vertical = side == BC_SIDE_LEFT || side == BC_SIDE_RIGHT;
    const int n = is_vertical ? this->rows : this->cols;

    this->boundaries[side] =
        std::vector<BoundaryElement<T>>(n, BoundaryElement<T>(value));
  }

  /**
   * @brief Specifically sets the boundary element of the specified side
   *        defined by the index to the boundary condition closed.
   *
   * @param side Side in which an element is to be defined as closed.
   * @param index Index of the boundary element on the corresponding
   *              boundary side. Must index an element of the corresponding
   * side.
   */
  void setBoundaryElemenClosed(BC_SIDE side, int index) {
    // tests whether the index really points to an element of the boundary side.
    tug_assert(boundaries[side].size() > index && index >= 0,
               "Index is selected either too large or too small.");

    this->boundaries[side][index].setType(BC_TYPE_CLOSED);
  }

  /**
   * @brief Specifically sets the boundary element of the specified side
   *        defined by the index to the boundary condition constant with the
            given concentration value.
   *
   * @param side Side in which an element is to be defined as constant.
   * @param index Index of the boundary element on the corresponding
   *              boundary side. Must index an element of the corresponding
   side.
   * @param value Concentration value to which the boundary element should be
   set.
   */
  void setBoundaryElementConstant(BC_SIDE side, int index, double value) {
    // tests whether the index really points to an element of the boundary side.
    tug_assert(boundaries[side].size() > index && index >= 0,
               "Index is selected either too large or too small.");
    this->boundaries[side][index].setType(BC_TYPE_CONSTANT);
    this->boundaries[side][index].setValue(value);
  }

  /**
   * @brief Returns the boundary condition of a specified side as a vector
   *        of BoundarsElement objects.
   *
   * @param side Boundary side from which the boundary conditions are to be
   * returned.
   * @return Contains the boundary conditions as
   * BoundaryElement<T> objects.
   */
  const std::vector<BoundaryElement<T>> &getBoundarySide(BC_SIDE side) const {
    tug_assert((this->dim > 1) ||
                   ((side == BC_SIDE_LEFT) || (side == BC_SIDE_RIGHT)),
               "For the "
               "one-dimensional case, only the BC_SIDE_LEFT and BC_SIDE_RIGHT "
               "borders exist.");
    return this->boundaries[side];
  }

  /**
   * @brief Get thes Boundary Side Values as a vector. Value is -1 in case some
   specific boundary is closed.
   *
   * @param side Boundary side for which the values are to be returned.
   * @return Vector with values as T.
   */
  Eigen::VectorX<T> getBoundarySideValues(BC_SIDE side) const {
    const std::size_t length = boundaries[side].size();
    Eigen::VectorX<T> values(length);

    for (int i = 0; i < length; i++) {
      if (getBoundaryElementType(side, i) == BC_TYPE_CLOSED) {
        values(i) = -1;
        continue;
      }
      values(i) = getBoundaryElementValue(side, i);
    }

    return values;
  }

  /**
   * @brief Returns the boundary condition of a specified element on a given
   * side.
   *
   * @param side Boundary side in which the boundary condition is located.
   * @param index Index of the boundary element on the corresponding
   *              boundary side. Must index an element of the corresponding
   * side.
   * @return Boundary condition as a BoundaryElement<T>
   * object.
   */
  BoundaryElement<T> getBoundaryElement(BC_SIDE side, int index) const {
    tug_assert(boundaries[side].size() > index && index >= 0,
               "Index is selected either too large or too small.");
    return this->boundaries[side][index];
  }

  /**
   * @brief Returns the type of a boundary condition, i.e. either BC_TYPE_CLOSED
   or BC_TYPE_CONSTANT.
   *
   * @param side Boundary side in which the boundary condition type is located.
   * @param index Index of the boundary element on the corresponding
   *              boundary side. Must index an element of the corresponding
   side.
   * @return Boundary Type of the corresponding boundary condition.
   */
  BC_TYPE getBoundaryElementType(BC_SIDE side, int index) const {
    tug_assert(boundaries[side].size() > index && index >= 0,
               "Index is selected either too large or too small.");
    return this->boundaries[side][index].getType();
  }

  /**
   * @brief Returns the concentration value of a corresponding
   *        BoundaryElement<T> object if it is a constant boundary condition.
   *
   * @param side Boundary side in which the boundary condition value is
   *             located.
   * @param index Index of the boundary element on the corresponding
   *              boundary side. Must index an element of the corresponding
   *              side.
   * @return Concentration of the corresponding BoundaryElement<T>
   * object.
   */
  T getBoundaryElementValue(BC_SIDE side, int index) const {
    tug_assert(boundaries[side].size() > index && index >= 0,
               "Index is selected either too large or too small.");
    tug_assert(
        boundaries[side][index].getType() == BC_TYPE_CONSTANT,
        "A value can only be output if it is a constant boundary condition.");

    return this->boundaries[side][index].getValue();
  }

  /**
   * @brief Serializes the boundary conditions into a vector of BoundaryElement
   * objects.
   *
   * @return Vector with BoundaryElement objects.
   */
  std::vector<BoundaryElement<T>> serialize() const {
    std::vector<BoundaryElement<T>> serialized;
    for (std::size_t side = 0; side < boundaries.size(); side++) {
      for (std::size_t index = 0; index < boundaries[side].size(); index++) {
        serialized.push_back(boundaries[side][index]);
      }
    }
    return serialized;
  }

  /**
   * @brief Deserializes the boundary conditions from a vector of
   * BoundaryElement objects.
   *
   * @param serialized Vector with BoundaryElement objects.
   */
  void deserialize(const std::vector<BoundaryElement<T>> &serialized) {
    std::size_t index = 0;
    for (std::size_t side = 0; side < boundaries.size(); side++) {
      for (std::size_t i = 0; i < boundaries[side].size(); i++) {
        boundaries[side][i] = serialized[index];
        index++;
      }
    }
  }

  /**
   *
   * @param index Index of the inner constant boundary condition
   * @param value Value of the inner constant boundary condition
   */
  void setInnerBoundary(std::uint32_t index, T value) {
    tug_assert(this->dim == 1, "This function is only available for 1D grids.");
    tug_assert(index < this->cols, "Index is out of bounds.");

    this->inner_boundary[std::make_pair(0, index)] = value;
  }

  /**
   * @brief Set inner constant boundary condition in 2D case.
   *
   * @param row Row index of the inner constant boundary condition
   * @param col Column index of the inner constant boundary condition
   * @param value Value of the inner constant boundary condition
   */
  void setInnerBoundary(std::uint32_t row, std::uint32_t col, T value) {
    tug_assert(this->dim == 2, "This function is only available for 2D grids.");
    tug_assert(row < this->rows && col < this->cols, "Index is out of bounds.");

    this->inner_boundary[std::make_pair(row, col)] = value;
  }

  /**
   * @brief Get inner constant boundary condition in 1D case.
   *
   * @param index Index of the inner constant boundary condition
   * @return std::pair<bool, T> Pair of boolean (whether constant boundary was
   * set or not) and value of the inner constant boundary condition
   */
  std::pair<bool, T> getInnerBoundary(std::uint32_t index) const {
    tug_assert(this->dim == 1, "This function is only available for 1D grids.");
    tug_assert(index < this->cols, "Index is out of bounds.");

    auto it = this->inner_boundary.find(std::make_pair(0, index));
    if (it == this->inner_boundary.end()) {
      return std::make_pair(false, -1);
    }
    return std::make_pair(true, it->second);
  }

  /**
   * @brief Get inner constant boundary condition in 2D case.
   *
   * @param row Row index of the inner constant boundary condition
   * @param col Column index of the inner constant boundary condition
   * @return std::pair<bool, T> Pair of boolean (whether constant boundary was
   * set or not) and value of the inner constant boundary condition
   */
  std::pair<bool, T> getInnerBoundary(std::uint32_t row,
                                      std::uint32_t col) const {
    tug_assert(this->dim == 2, "This function is only available for 2D grids.");
    tug_assert(row < this->rows && col < this->cols, "Index is out of bounds.");

    auto it = this->inner_boundary.find(std::make_pair(row, col));
    if (it == this->inner_boundary.end()) {
      return std::make_pair(false, -1);
    }
    return std::make_pair(true, it->second);
  }

  /**
   * @brief Get inner constant boundary conditions of a row as a vector. Can be
   * used for 1D grids (row == 0) or 2D grids.
   *
   * @param row Index of the row for which the inner boundary conditions are to
   * be returned.
   * @return std::vector<std::pair<bool, T>> Vector of pairs of boolean (whether
   * constant boundary was set or not) and value of the inner constant boundary
   * condition
   */
  std::vector<std::pair<bool, T>> getInnerBoundaryRow(std::uint32_t row) const {
    tug_assert(row < this->rows, "Index is out of bounds.");

    if (inner_boundary.empty()) {
      return std::vector<std::pair<bool, T>>(this->cols,
                                             std::make_pair(false, -1));
    }

    std::vector<std::pair<bool, T>> row_values;
    for (std::uint32_t col = 0; col < this->cols; col++) {
      row_values.push_back(getInnerBoundary(row, col));
    }

    return row_values;
  }

  /**
   * @brief Get inner constant boundary conditions of a column as a vector. Can
   * only be used for 2D grids.
   *
   * @param col Index of the column for which the inner boundary conditions are
   * to be returned.
   * @return std::vector<std::pair<bool, T>> Vector of pairs of boolean (whether
   * constant boundary was set or not) and value of the inner constant boundary
   * condition
   */
  std::vector<std::pair<bool, T>> getInnerBoundaryCol(std::uint32_t col) const {
    tug_assert(this->dim == 2, "This function is only available for 2D grids.");
    tug_assert(col < this->cols, "Index is out of bounds.");

    if (inner_boundary.empty()) {
      return std::vector<std::pair<bool, T>>(this->rows,
                                             std::make_pair(false, -1));
    }

    std::vector<std::pair<bool, T>> col_values;
    for (std::uint32_t row = 0; row < this->rows; row++) {
      col_values.push_back(getInnerBoundary(row, col));
    }

    return col_values;
  }

  /**
   * @brief Get inner constant boundary conditions as a map. Can be read by
   * setInnerBoundaries.
   *
   * @return Map of inner constant boundary conditions
   */
  std::map<std::pair<std::uint32_t, std::uint32_t>, T>
  getInnerBoundaries() const {
    return this->inner_boundary;
  }

  /**
   * @brief Set inner constant boundary conditions as a map. Can be obtained by
   * getInnerBoundaries.
   *
   * @param inner_boundary Map of inner constant boundary conditions
   */
  void
  setInnerBoundaries(const std::map<std::pair<std::uint32_t, std::uint32_t>, T>
                         &inner_boundary) {
    this->inner_boundary = inner_boundary;
  }

private:
  const std::uint8_t dim;
  const std::uint32_t cols;
  const std::uint32_t rows;

  std::vector<std::vector<BoundaryElement<T>>>
      boundaries; // Vector with Boundary Element information

  // Inner boundary conditions for 1D and 2D grids identified by (row,col)
  std::map<std::pair<std::uint32_t, std::uint32_t>, T> inner_boundary;
};
} // namespace tug
#endif // BOUNDARY_H_
