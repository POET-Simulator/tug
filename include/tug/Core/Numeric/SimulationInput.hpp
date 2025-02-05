#pragma once

#include <cstdint>
#include <tug/Boundary.hpp>
#include <tug/Core/Matrix.hpp>

namespace tug {

template <typename T> struct SimulationInput {
  RowMajMatMap<T> &concentrations;
  const RowMajMat<T> &alphaX;
  const RowMajMat<T> &alphaY;
  const Boundary<T> boundaries;

  const std::uint8_t dim;
  const T timestep;
  const std::size_t rowMax;
  const std::size_t colMax;
  const T deltaRow;
  const T deltaCol;
};
} // namespace tug
