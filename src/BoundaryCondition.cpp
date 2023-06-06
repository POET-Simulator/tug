#include <algorithm>
#include <tug/BoundaryCondition.hpp>
#include <vector>

#include "TugUtils.hpp"

constexpr uint8_t DIM_1D = 2;
constexpr uint8_t DIM_2D = 4;

tug::bc::BoundaryCondition::BoundaryCondition(int x) {
  this->bc_internal.resize(DIM_1D, {0, 0});
  this->dim = 1;
  // this value is actually unused
  this->maxsize = 1;

  this->sizes[X_DIM] = x;
  this->sizes[Y_DIM] = 1;

  this->maxindex = x - 1;
}

tug::bc::BoundaryCondition::BoundaryCondition(int x, int y) {
  this->maxsize = (x >= y ? x : y);
  this->bc_internal.resize(DIM_2D * maxsize, {0, 0});
  this->dim = 2;

  this->sizes[X_DIM] = x;
  this->sizes[Y_DIM] = y;

  this->maxindex = (x * y) - 1;
}

void tug::bc::BoundaryCondition::setSide(
    uint8_t side, tug::bc::boundary_condition &input_bc) {
  // QUESTION: why cant the BC be changed for dim = 1?
  if (this->dim == 1) {
    throw_invalid_argument("setSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw_out_of_range("Invalid range for 2D grid");
  }

  uint32_t size =
      (side == tug::bc::BC_SIDE_LEFT || side == tug::bc::BC_SIDE_RIGHT
           ? this->sizes[Y_DIM]
           : this->sizes[X_DIM]);

  for (uint32_t i = 0; i < size; i++) {
    this->bc_internal[side * maxsize + i] = input_bc;
  }
}

void tug::bc::BoundaryCondition::setSide(
    uint8_t side, std::vector<tug::bc::boundary_condition> &input_bc) {
  if (this->dim == 1) {
    throw_invalid_argument("setSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw_out_of_range("Invalid range for 2D grid");
  }

  uint32_t size =
      (side == tug::bc::BC_SIDE_LEFT || side == tug::bc::BC_SIDE_RIGHT
           ? this->sizes[Y_DIM]
           : this->sizes[X_DIM]);

  if (input_bc.size() > size) {
    throw_out_of_range("Input vector is greater than maximum excpected value");
  }

  for (int i = 0; i < size; i++) {
    bc_internal[this->maxsize * side + i] = input_bc[i];
  }
}

auto tug::bc::BoundaryCondition::getSide(uint8_t side)
    -> std::vector<tug::bc::boundary_condition> {
  if (this->dim == 1) {
    throw_invalid_argument("getSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw_out_of_range("Invalid range for 2D grid");
  }

  uint32_t size =
      (side == tug::bc::BC_SIDE_LEFT || side == tug::bc::BC_SIDE_RIGHT
           ? this->sizes[Y_DIM]
           : this->sizes[X_DIM]);

  std::vector<tug::bc::boundary_condition> out(size);

  for (int i = 0; i < size; i++) {
    out[i] = this->bc_internal[this->maxsize * side + i];
  }

  return out;
}

auto tug::bc::BoundaryCondition::row_boundary(uint32_t i) const
    -> tug::bc::bc_tuple {
  if (i >= this->sizes[Y_DIM]) {
    throw_out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_LEFT * this->maxsize + i],
          this->bc_internal[BC_SIDE_RIGHT * this->maxsize + i]};
}

auto tug::bc::BoundaryCondition::col_boundary(uint32_t i) const
    -> tug::bc::bc_tuple {
  if (this->dim == 1) {
    throw_invalid_argument("Access of column requires at least 2D grid");
  }
  if (i >= this->sizes[X_DIM]) {
    throw_out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_TOP * this->maxsize + i],
          this->bc_internal[BC_SIDE_BOTTOM * this->maxsize + i]};
}

auto tug::bc::BoundaryCondition::getInnerRow(uint32_t i) const -> bc_vec {
  if (i >= this->sizes[Y_DIM]) {
    throw_out_of_range("Index is out of range");
  }

  bc_vec row(this->sizes[X_DIM], {tug::bc::BC_UNSET, 0});

  if (this->inner_cells.empty()) {
    return row;
  }

  uint32_t index_min = i * this->sizes[X_DIM];
  uint32_t index_max = ((i + 1) * this->sizes[X_DIM]) - 1;

  for (auto const &cell : this->inner_cells) {
    if (cell.first < index_min) {
      continue;
    }
    if (cell.first > index_max) {
      break;
    }

    row[cell.first - index_min] = cell.second;
  }

  return row;
}

auto tug::bc::BoundaryCondition::getInnerCol(uint32_t i) const -> bc_vec {
  if (this->dim != 2) {
    throw_invalid_argument("getInnerCol is only applicable for 2D grids");
  }
  if (i >= this->sizes[X_DIM]) {
    throw_out_of_range("Index is out of range");
  }

  bc_vec col(this->sizes[Y_DIM], {tug::bc::BC_UNSET, 0});

  if (this->inner_cells.empty()) {
    return col;
  }

  for (auto const &cell : this->inner_cells) {
    if (cell.first % this->sizes[X_DIM] == i) {
      col[cell.first / this->sizes[X_DIM]] = cell.second;
    }
  }

  return col;
}

void tug::bc::BoundaryCondition::setInnerBC(boundary_condition bc, int x,
                                            int y = 0) {
  if (x >= this->sizes[X_DIM] || y >= this->sizes[Y_DIM]) {
    throw_out_of_range("One input parameter is out of range");
  }
  uint32_t index = y * this->sizes[X_DIM] + x;
  auto it = this->inner_cells.find(index);

  if (it != this->inner_cells.end()) {
    it->second = bc;
    return;
  }

  this->inner_cells.insert({index, bc});
}

void tug::bc::BoundaryCondition::unsetInnerBC(int x, int y) {
  uint32_t index = y * this->sizes[X_DIM] + x;
  this->inner_cells.erase(index);
}

auto tug::bc::BoundaryCondition::getInnerBC(int x, int y = 0)
    -> boundary_condition {
  if (x >= this->sizes[X_DIM] || y >= this->sizes[Y_DIM]) {
    throw_out_of_range("One input parameter is out of range");
  }

  uint32_t index = y * this->sizes[X_DIM] + x;

  auto it = this->inner_cells.find(index);

  if (it != this->inner_cells.end()) {
    return it->second;
  }

  return {tug::bc::BC_UNSET, 0};
}
