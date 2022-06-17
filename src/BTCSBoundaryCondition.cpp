#include <algorithm>
#include <bits/stdint-uintn.h>
#include <diffusion/BTCSBoundaryCondition.hpp>
#include <stdexcept>
#include <vector>

constexpr uint8_t DIM_1D = 2;
constexpr uint8_t DIM_2D = 4;

Diffusion::BTCSBoundaryCondition::BTCSBoundaryCondition() {
  this->bc_internal.resize(DIM_1D, {0, 0});
  this->dim = 1;
  // this value is actually unused
  this->maxsize = 1;
}

Diffusion::BTCSBoundaryCondition::BTCSBoundaryCondition(int x, int y) {
  this->maxsize = (x >= y ? x : y);
  this->bc_internal.resize(DIM_2D * maxsize, {0, 0});
  this->dim = 2;
}

void Diffusion::BTCSBoundaryCondition::setSide(
    uint8_t side, Diffusion::boundary_condition &input_bc) {
  if (this->dim == 1) {
    throw std::invalid_argument("setSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw std::out_of_range("Invalid range for 2D grid");
  }

  for (int i = 0; i < maxsize; i++) {
    this->bc_internal[side * maxsize + i] = input_bc;
  }
}

void Diffusion::BTCSBoundaryCondition::setSide(
    uint8_t side, std::vector<Diffusion::boundary_condition> &input_bc) {
  if (this->dim == 1) {
    throw std::invalid_argument("setSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw std::out_of_range("Invalid range for 2D grid");
  }
  if (input_bc.size() > this->maxsize) {
    throw std::out_of_range(
        "Input vector is greater than maximum excpected value");
  }

  for (int i = 0; i < input_bc.size(); i++) {
    bc_internal[this->maxsize * side + i] = input_bc[i];
  }
}

auto Diffusion::BTCSBoundaryCondition::getSide(uint8_t side)
    -> std::vector<Diffusion::boundary_condition> {
  if (this->dim == 1) {
    throw std::invalid_argument("getSide requires at least a 2D grid");
  }
  if (side > 3) {
    throw std::out_of_range("Invalid range for 2D grid");
  }

  std::vector<Diffusion::boundary_condition> out(this->maxsize);

  for (int i = 0; i < this->maxsize; i++) {
    out[i] = this->bc_internal[this->maxsize * side + i];
  }

  return out;
}

auto Diffusion::BTCSBoundaryCondition::col(uint32_t i) const
    -> Diffusion::bc_tuple {
  if (this->dim == 1) {
    throw std::invalid_argument("Access of column requires at least 2D grid");
  }
  if (i >= this->maxsize) {
    throw std::out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_TOP * this->maxsize + i],
          this->bc_internal[BC_SIDE_BOTTOM * this->maxsize + i]};
}

auto Diffusion::BTCSBoundaryCondition::row(uint32_t i) const
    -> Diffusion::bc_tuple {
  if (i >= this->maxsize) {
    throw std::out_of_range("Index out of range");
  }

  return {this->bc_internal[BC_SIDE_LEFT * this->maxsize + i],
          this->bc_internal[BC_SIDE_RIGHT * this->maxsize + i]};
}
