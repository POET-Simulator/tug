/**
 * @file Boundary.hpp
 * @brief API of Boundary class, that holds all information for each boundary
 * condition at the edges of the diffusion grid.
 *
 */
#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "Grid.hpp"
#include <cstddef>
#include <cstdint>
#include <exception>

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
    if (value < 0) {
      throw std::invalid_argument("No negative concentration allowed.");
    }
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
    if (this->dim == 1) {
      if ((side == BC_SIDE_BOTTOM) || (side == BC_SIDE_TOP)) {
        throw std::invalid_argument(
            "For the one-dimensional case, only the BC_SIDE_LEFT and "
            "BC_SIDE_RIGHT borders exist.");
      }
    }

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
    if (this->dim == 1) {
      if ((side == BC_SIDE_BOTTOM) || (side == BC_SIDE_TOP)) {
        throw std::invalid_argument(
            "For the one-dimensional case, only the BC_SIDE_LEFT and "
            "BC_SIDE_RIGHT borders exist.");
      }
    }

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
    if ((boundaries[side].size() < index) || index < 0) {
      throw std::invalid_argument(
          "Index is selected either too large or too small.");
    }
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
    if ((boundaries[side].size() < index) || index < 0) {
      throw std::invalid_argument(
          "Index is selected either too large or too small.");
    }
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
    if (this->dim == 1) {
      if ((side == BC_SIDE_BOTTOM) || (side == BC_SIDE_TOP)) {
        throw std::invalid_argument(
            "For the one-dimensional trap, only the BC_SIDE_LEFT and "
            "BC_SIDE_RIGHT borders exist.");
      }
    }
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
    if ((boundaries[side].size() < index) || index < 0) {
      throw std::invalid_argument(
          "Index is selected either too large or too small.");
    }
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
    if ((boundaries[side].size() < index) || index < 0) {
      throw std::invalid_argument(
          "Index is selected either too large or too small.");
    }
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
    if ((boundaries[side].size() < index) || index < 0) {
      throw std::invalid_argument(
          "Index is selected either too large or too small.");
    }
    if (boundaries[side][index].getType() != BC_TYPE_CONSTANT) {
      throw std::invalid_argument(
          "A value can only be output if it is a constant boundary condition.");
    }
    return this->boundaries[side][index].getValue();
  }

private:
  const std::uint8_t dim;
  const std::uint32_t cols;
  const std::uint32_t rows;

  std::vector<std::vector<BoundaryElement<T>>>
      boundaries; // Vector with Boundary Element information
};
} // namespace tug
#endif // BOUNDARY_H_
