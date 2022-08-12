#include <algorithm>
#include <bits/stdint-uintn.h>
#include <vector>

#include "grid/BTCSBoundaryCondition.hpp"
#include "diffusion/BTCSUtils.hpp"

constexpr uint8_t DIM_1D = 2;
constexpr uint8_t DIM_2D = 4;

Diffusion::BTCSBoundaryCondition::BTCSBoundaryCondition() {
  this->bc_internal.resize(DIM_1D, {0, 0});
  this->dim = 1;
  // this value is actually unused
  this->maxsize = 1;

  this->sizes[0] = 1;
  this->sizes[1] = 0;
}

Diffusion::BTCSBoundaryCondition::BTCSBoundaryCondition(int x, int y) {
  this->maxsize = (x >= y ? x : y);
  this->bc_internal.resize(DIM_2D * maxsize, {0, 0});
  this->dim = 2;

  this->sizes[0] = x;
  this->sizes[1] = y;
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
           ? this->sizes[0]
           : this->sizes[1]);

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
           ? this->sizes[0]
           : this->sizes[1]);

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
           ? this->sizes[0]
           : this->sizes[1]);

  std::vector<Diffusion::boundary_condition> out(size);

  for (int i = 0; i < size; i++) {
    out[i] = this->bc_internal[this->maxsize * side + i];
  }

  return out;
}

auto Diffusion::BTCSBoundaryCondition::col(uint32_t i) const
    -> Diffusion::bc_tuple {
  if (this->dim == 1) {
    throw_invalid_argument("Access of column requires at least 2D grid");
  }
  if (i >= this->sizes[1]) {
    throw_out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_TOP * this->maxsize + i],
          this->bc_internal[BC_SIDE_BOTTOM * this->maxsize + i]};
}

auto Diffusion::BTCSBoundaryCondition::row(uint32_t i) const
    -> Diffusion::bc_tuple {
  if (i >= this->sizes[0]) {
    throw_out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_LEFT * this->maxsize + i],
          this->bc_internal[BC_SIDE_RIGHT * this->maxsize + i]};
}
