#ifndef BOUNDARYCONDITION_H_
#define BOUNDARYCONDITION_H_

#include <array>
#include <bits/stdint-uintn.h>
#include <stdexcept>
#include <stdint.h>
#include <vector>

typedef uint8_t bctype;

namespace tug {
namespace boundary_condition {

enum {
  BC_TYPE_CLOSED,   /**< Defines a closed/Neumann boundary condition. */
  BC_TYPE_FLUX,     /**< Defines a flux/Cauchy boundary condition. */
  BC_TYPE_CONSTANT, /**< Defines a constant/Dirichlet boundary condition. */
  BC_UNSET          /**< Indicates undefined boundary condition*/
};

enum {
  BC_SIDE_LEFT,  /**< Defines boundary conditions for the left side of the grid.
                  */
  BC_SIDE_RIGHT, /**< Defines boundary conditions for the right side of the
                    grid. */
  BC_SIDE_TOP,   /**< Defines boundary conditions for the top of the grid. */
  BC_SIDE_BOTTOM, /**< Defines boundary conditions for the bottom of the grid.
                   */
  BC_INNER
};

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
 * Represents both boundary conditions of a row/column.
 */
typedef std::array<boundary_condition, 2> bc_tuple;
typedef std::vector<boundary_condition> bc_vec;

/**
 * Class to define the boundary condition of a grid.
 */
class BTCSBoundaryCondition {
public:
  /**
   *  Creates a new instance with two elements. Used when defining boundary
   *  conditions of 1D grids.
   *
   *  \param x Number of grid cells in x-direction
   */
  BTCSBoundaryCondition(int x);

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
  auto row_boundary(uint32_t i) const -> bc_tuple;

  /**
   * Get both boundary conditions of a given column (top and bottom).
   *
   * \param i Index of column
   *
   * \returns Top and bottom boundary values as an array (defined as data
   * type bc_tuple).
   */
  auto col_boundary(uint32_t i) const -> bc_tuple;

  /**
   * Get a row of field and its inner boundary conditions.
   *
   * \param i Index of the row starting at 0.
   *
   * \returns Row of the inner boundary conditions of the field.
   */
  auto getInnerRow(uint32_t i) const -> bc_vec;

  /**
   * Get a column of field and its inner boundary conditions.
   *
   * \param i Index of the column starting at 0.
   *
   * \returns Column of the inner boundary conditions of the field.
   */
  auto getInnerCol(uint32_t i) const -> bc_vec;

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

  bc_vec special_cells;

  uint8_t dim;

  uint32_t sizes[2];
  uint32_t maxsize;
  uint32_t maxindex;

  enum { X_DIM, Y_DIM };

public:
  /**
   * Returns the left/right boundary condition for 1D grid.
   *
   * \param side Side of the boundary condition to get.
   *
   * \returns Boundary condition
   */
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

  /**
   * Returns the boundary condition of a given side for 2D grids.
   *
   * \param side Side of the boundary condition to get.
   * \param i Index of the boundary condition.
   *
   * \returns Boundary condition
   */
  boundary_condition operator()(uint8_t side, uint32_t i) const {
    if (side == BC_INNER) {
      if (i > maxindex) {
        throw std::out_of_range("Index exceeds grid cell numbers");
      }
      return special_cells[i];
    }
    if (dim != 2) {
      throw std::invalid_argument(
          "Only 2D grids support 2 parameters in operator");
    }
    if (side > 3) {
      throw std::out_of_range("2D index out of range");
    }
    return bc_internal[side * maxsize + i];
  }

  /**
   * Returns the left/right boundary condition for 1D grid.
   *
   * \param side Side of the boundary condition to get.
   *
   * \returns Boundary condition
   */
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

  /**
   * Returns the boundary condition of a given side for 2D grids.
   *
   * \param side Side of the boundary condition to get.
   * \param i Index of the boundary condition.
   *
   * \returns Boundary condition
   */
  boundary_condition &operator()(uint8_t side, uint32_t i) {
    if (side == BC_INNER) {
      if (i > maxindex) {
        throw std::out_of_range("Index exceeds grid cell numbers");
      }
      return special_cells[i];
    }
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

} // namespace boundary_condition
} // namespace tug
#endif // BOUNDARYCONDITION_H_
