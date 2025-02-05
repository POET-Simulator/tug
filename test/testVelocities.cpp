#include "tug/Boundary.hpp"
#include "tug/Core/Matrix.hpp"
#include <cstddef>
#include <tug/Advection/Velocities.hpp>

#include <gtest/gtest.h>

#define VELOCITIES_TEST(x) TEST(Velocities, x)

VELOCITIES_TEST(SteadyStateCenter) {
  // Arrange
  constexpr std::size_t rows = 25;
  constexpr std::size_t cols = 25;

  constexpr std::size_t center_row = rows / 2;
  constexpr std::size_t center_col = cols / 2;

  constexpr double K = 1E-2;

  tug::RowMajMat<double> hydHeads =
      tug::RowMajMat<double>::Constant(rows, cols, 1);

  tug::RowMajMat<double> permKX =
      tug::RowMajMat<double>::Constant(rows, cols, K);
  tug::RowMajMat<double> permKY =
      tug::RowMajMat<double>::Constant(rows, cols, K);

  tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                  tug::HYDRAULIC_RESOLVE::EXPLICIT>
      velo(hydHeads);

  velo.setDomain(100, 100);
  velo.setAlphaX(permKX);
  velo.setAlphaY(permKY);

  tug::Boundary<double> &bcH = velo.getBoundaryConditions();
  bcH.setBoundarySideConstant(tug::BC_SIDE_LEFT, 1);
  bcH.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 1);
  bcH.setBoundarySideConstant(tug::BC_SIDE_TOP, 1);
  bcH.setBoundarySideConstant(tug::BC_SIDE_BOTTOM, 1);

  bcH.setInnerBoundary(center_row, center_col, 10);

  velo.run();

  const auto &velocitiesX = velo.getVelocitiesX();
  const auto &velocitiesY = velo.getVelocitiesY();

  // Assert

  // check velocities in x-direction
  for (std::size_t i_rows = 0; i_rows < rows; i_rows++) {
    for (std::size_t i_cols = 0; i_cols < cols + 1; i_cols++) {
      if (i_rows <= center_row && i_cols <= center_col) {
        EXPECT_LE(velocitiesX(i_rows, i_cols), 0);
      } else if (i_rows > center_row && i_cols > center_col) {
        EXPECT_GE(velocitiesX(i_rows, i_cols), 0);
      } else if (i_rows <= center_row && i_cols > center_col) {
        EXPECT_GE(velocitiesX(i_rows, i_cols), 0);
      } else if (i_rows > center_row && i_cols <= center_col) {
        EXPECT_LE(velocitiesX(i_rows, i_cols), 0);
      } else {
        FAIL() << "Uncovered case";
      }
    }
  }

  // check velocities in y-direction
  for (std::size_t i_rows = 0; i_rows < rows + 1; i_rows++) {
    for (std::size_t i_cols = 0; i_cols < cols; i_cols++) {
      if (i_rows <= center_row && i_cols <= center_col) {
        EXPECT_LE(velocitiesY(i_rows, i_cols), 0);
      } else if (i_rows > center_row && i_cols > center_col) {
        EXPECT_GE(velocitiesY(i_rows, i_cols), 0);
      } else if (i_rows <= center_row && i_cols > center_col) {
        EXPECT_LE(velocitiesY(i_rows, i_cols), 0);
      } else if (i_rows > center_row && i_cols <= center_col) {
        EXPECT_GE(velocitiesY(i_rows, i_cols), 0);
      } else {
        FAIL() << "Uncovered case";
      }
    }
  }
}
