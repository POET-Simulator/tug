#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include <array>
#include <bits/stdint-uintn.h>
#include <stdexcept>
#include <stdint.h>
#include <vector>

typedef uint8_t bctype;

namespace Diffusion {

/**
 * Defines a closed/Neumann boundary condition.
 */
static const bctype BC_TYPE_CLOSED = 0;

/**
 * Defines a flux/Cauchy boundary condition.
 */
static const bctype BC_TYPE_FLUX = 1;

/**
 * Defines a constant/Dirichlet boundary condition.
 */
static const bctype BC_TYPE_CONSTANT = 2;

/**
 * Defines boundary conditions for the left side of the grid.
 */
static const uint8_t BC_SIDE_LEFT = 0;

/**
 * Defines boundary conditions for the right side of the grid.
 */
static const uint8_t BC_SIDE_RIGHT = 1;

/**
 * Defines boundary conditions for the top of the grid.
 */
static const uint8_t BC_SIDE_TOP = 2;

/**
 * Defines boundary conditions for the bottom of the grid.
 */
static const uint8_t BC_SIDE_BOTTOM = 3;

/**
 * Defines the boundary condition type and according value.
 */
typedef struct boundary_condition_s {
  bctype type;  /**< Type of the boundary condition */
  double value; /**< Value of the boundary condition. Either a concrete value of
                   concentration for BC_TYPE_CONSTANT or gradient when type is
                   BC_TYPE_FLUX. Unused if BC_TYPE_CLOSED.*/
} boundary_condition;

/**
 * Internal use only.
 */
typedef std::array<boundary_condition, 2> bc_tuple;

/**
 * Class to define the boundary condition of a grid.
 */
class BTCSBoundaryCondition {
public:
  /**
   *  Creates a new instance with two elements. Used when defining boundary
   *  conditions of 1D grids.
   */
  BTCSBoundaryCondition();

  /**
   * Creates a new instance with 4 * max(x,y) elements. Used to describe the
   * boundary conditions for 2D grids. NOTE: On non-squared grids there are more
   * elements than needed and no exception is thrown if some index exceeding
   * grid limits.
   *
   * \param x Number of grid cells in x-direction
   * \param y Number of grid cells in y-direction
   *
   */
  BTCSBoundaryCondition(int x, int y);

  /**
   * Sets the boundary condition for a specific side of the grid.
   *
   * \param side Side for which the given boundary condition should be set.
   * \param input_bc Instance of struct boundary_condition with desired boundary
   * condition.
   *
   * \throws std::invalid_argument Indicates wrong dimensions of the grid.
   * \throws std::out_of_range Indicates a out of range value for side.
   */
  void setSide(uint8_t side, boundary_condition &input_bc);

  /**
   * Sets the boundary condition for a specific side of the grid.
   *
   * \param side Side for which the given boundary condition should be set.
   * \param input_bc Vector of boundary conditions for specific side.
   *
   * \throws std::invalid_argument Indicates wrong dimensions of the grid.
   * \throws std::out_of_range Indicates a out of range value for side or
   * invalid size of input vector.
   */
  void setSide(uint8_t side, std::vector<boundary_condition> &input_bc);

  /**
   * Returns a vector of boundary conditions of given side. Can be used to set
   * custom boundary conditions and set back via setSide() with vector input.
   *
   * \param side Side which boundary conditions should be returned
   *
   * \returns Vector of boundary conditions
   *
   * \throws std::invalid_argument If given dimension is less or equal to 1.
   * \throws std::out_of_range Indicates a out of range value for side.
   */
  auto getSide(uint8_t side) -> std::vector<boundary_condition>;

  /**
   * Get both boundary conditions of a given row (left and right).
   *
   * \param i Index of row
   *
   * \returns Left and right boundary values as an array (defined as data
   * type bc_tuple).
   */
  auto row(uint32_t i) const -> bc_tuple;

  /**
   * Get both boundary conditions of a given column (top and bottom).
   *
   * \param i Index of column
   *
   * \returns Top and bottom boundary values as an array (defined as data
   * type bc_tuple).
   */
  auto col(uint32_t i) const -> bc_tuple;

  /**
   * Create an instance of boundary_condition data type. Can be seen as a helper
   * function.
   *
   * \param type Type of the boundary condition.
   * \param value According value of condition.
   *
   * \returns Instance of boundary_condition
   */
  static boundary_condition returnBoundaryCondition(bctype type, double value) {
    return {type, value};
  }

private:
  std::vector<boundary_condition> bc_internal;

  uint8_t dim;

  uint32_t sizes[2];
  uint32_t maxsize;

public:
  boundary_condition operator()(uint8_t side) const {
    if (dim != 1) {
      throw std::invalid_argument(
          "Only 1D grid support 1 parameter in operator");
    }
    if (side > 1) {
      throw std::out_of_range("1D index out of range");
    }
    return bc_internal[side];
  }

  boundary_condition operator()(uint8_t side, uint32_t i) const {
    if (dim != 2) {
      throw std::invalid_argument(
          "Only 2D grids support 2 parameters in operator");
    }
    if (side > 3) {
      throw std::out_of_range("2D index out of range");
    }
    return bc_internal[side * maxsize + i];
  }

  boundary_condition &operator()(uint8_t side) {
    if (dim != 1) {
      throw std::invalid_argument(
          "Only 1D grid support 1 parameter in operator");
    }
    if (side > 1) {
      throw std::out_of_range("1D index out of range");
    }
    return bc_internal[side];
  }

  boundary_condition &operator()(uint8_t side, uint32_t i) {
    if (dim != 2) {
      throw std::invalid_argument("Explicit setting of bc value with 2 "
                                  "parameters is only supported for 2D grids");
    }
    if (side > 3) {
      throw std::out_of_range("2D index out of range");
    }
    return bc_internal[side * maxsize + i];
  }
};

} // namespace Diffusion

#endif // BOUNDARYCONDITION_H_
