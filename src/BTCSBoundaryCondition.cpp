#include <algorithm>
#include <bits/stdint-uintn.h>
#include <vector>

#include "grid/BTCSBoundaryCondition.hpp"
#include "BTCSUtils.hpp"

constexpr uint8_t DIM_1D = 2;
constexpr uint8_t DIM_2D = 4;

Diffusion::BTCSBoundaryCondition::BTCSBoundaryCondition(int x) {
  this->bc_internal.resize(DIM_1D, {0, 0});
  this->special_cells.resize(x, {BC_UNSET, 0});
  this->dim = 1;
  // this value is actually unused
  this->maxsize = 1;

  this->sizes[X_DIM] = 1;
  this->sizes[Y_DIM] = x;

  this->maxindex = x - 1;
}

Diffusion::BTCSBoundaryCondition::BTCSBoundaryCondition(int x, int y) {
  this->maxsize = (x >= y ? x : y);
  this->bc_internal.resize(DIM_2D * maxsize, {0, 0});
  this->special_cells.resize(x * y, {BC_UNSET, 0});
  this->dim = 2;

  this->sizes[X_DIM] = x;
  this->sizes[Y_DIM] = y;

  this->maxindex = (x * y) - 1;
}

void Diffusion::BTCSBoundaryCondition::setSide(
    uint8_t side, Diffusion::boundary_condition &input_bc) {
  if (this->dim == 1) {
    throw_invalid_argument("setSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw_out_of_range("Invalid range for 2D grid");
  }

  uint32_t size =
      (side == Diffusion::BC_SIDE_LEFT || side == Diffusion::BC_SIDE_RIGHT
           ? this->sizes[X_DIM]
           : this->sizes[Y_DIM]);

  for (uint32_t i = 0; i < size; i++) {
    this->bc_internal[side * maxsize + i] = input_bc;
  }
}

void Diffusion::BTCSBoundaryCondition::setSide(
    uint8_t side, std::vector<Diffusion::boundary_condition> &input_bc) {
  if (this->dim == 1) {
    throw_invalid_argument("setSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw_out_of_range("Invalid range for 2D grid");
  }

  uint32_t size =
      (side == Diffusion::BC_SIDE_LEFT || side == Diffusion::BC_SIDE_RIGHT
           ? this->sizes[X_DIM]
           : this->sizes[Y_DIM]);

  if (input_bc.size() > size) {
    throw_out_of_range("Input vector is greater than maximum excpected value");
  }

  for (int i = 0; i < size; i++) {
    bc_internal[this->maxsize * side + i] = input_bc[i];
  }
}

auto Diffusion::BTCSBoundaryCondition::getSide(uint8_t side)
    -> std::vector<Diffusion::boundary_condition> {
  if (this->dim == 1) {
    throw_invalid_argument("getSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw_out_of_range("Invalid range for 2D grid");
  }

  uint32_t size =
      (side == Diffusion::BC_SIDE_LEFT || side == Diffusion::BC_SIDE_RIGHT
           ? this->sizes[X_DIM]
           : this->sizes[Y_DIM]);

  std::vector<Diffusion::boundary_condition> out(size);

  for (int i = 0; i < size; i++) {
    out[i] = this->bc_internal[this->maxsize * side + i];
  }

  return out;
}

auto Diffusion::BTCSBoundaryCondition::col_boundary(uint32_t i) const
    -> Diffusion::bc_tuple {
  if (this->dim == 1) {
    throw_invalid_argument("Access of column requires at least 2D grid");
  }
  if (i >= this->sizes[Y_DIM]) {
    throw_out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_TOP * this->maxsize + i],
          this->bc_internal[BC_SIDE_BOTTOM * this->maxsize + i]};
}

auto Diffusion::BTCSBoundaryCondition::row_boundary(uint32_t i) const
    -> Diffusion::bc_tuple {
  if (i >= this->sizes[X_DIM]) {
    throw_out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_LEFT * this->maxsize + i],
          this->bc_internal[BC_SIDE_RIGHT * this->maxsize + i]};
}

auto Diffusion::BTCSBoundaryCondition::getInnerRow(uint32_t i) const -> bc_vec {
  if (i >= this->sizes[X_DIM]) {
    throw_out_of_range("Index is out of range");
  }

  bc_vec::const_iterator start =
      this->special_cells.begin() + (i * this->sizes[Y_DIM]);
  bc_vec::const_iterator end =
      this->special_cells.begin() + ((i + 1) * this->sizes[Y_DIM]);

  bc_vec row(start, end);

  return row;
}

auto Diffusion::BTCSBoundaryCondition::getInnerCol(uint32_t i) const -> bc_vec {
  if (this->dim != 2) {
    throw_invalid_argument("getInnerCol is only applicable for 2D grids");
  } else if (i >= this->sizes[X_DIM]) {
    throw_out_of_range("Index is out of range");
  }

  bc_vec col;
  col.reserve(this->sizes[X_DIM]);

  for (int j = 0; j < this->sizes[X_DIM]; i += this->sizes[Y_DIM], j++) {
    col[j] = this->special_cells[i];
  }

  return col;
}
