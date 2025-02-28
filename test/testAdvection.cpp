#include <gtest/gtest.h>
#include <tug/tug.hpp>

#define ADVECTION_TEST(x) TEST(Advection, x)

ADVECTION_TEST(LeftToRight) {
  constexpr std::size_t rows = 21;
  constexpr std::size_t cols = 21;

  constexpr double K = 1E-2;
  constexpr double timestep = 5039.05;
  constexpr std::size_t iterations = 21;
  constexpr double porosity = 0.2;

  constexpr double epsilon = 1E-13;

  tug::RowMajMat<double> hydHeads =
      tug::RowMajMat<double>::Constant(rows, cols, 1);
  tug::RowMajMat<double> concentrations =
      tug::RowMajMat<double>::Constant(rows, cols, 0);

  tug::RowMajMat<double> permK =
      tug::RowMajMat<double>::Constant(rows, cols, K);

  tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                  tug::HYDRAULIC_RESOLVE::IMPLICIT>
      velocities(hydHeads);
  velocities.setDomain(100, 100);
  velocities.setPermKX(permK);
  velocities.setPermKY(permK);
  velocities.setEpsilon(1E-8);

  tug::Advection advection(concentrations, velocities);

  advection.setPorosity(tug::RowMajMat<double>::Constant(rows, cols, porosity));
  advection.setIterations(iterations);
  advection.setTimestep(timestep);

  tug::Boundary<double> &bcH = velocities.getBoundaryConditions();
  bcH.setBoundarySideConstant(tug::BC_SIDE_LEFT, 10);
  bcH.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 1);

  tug::Boundary<double> &bcC = advection.getBoundaryConditions();
  bcC.setBoundarySideConstant(tug::BC_SIDE_LEFT, 0.1);
  bcC.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 1);

  advection.run();

  // check if the concentration is transported from left to right
  for (std::size_t i_rows = 0; i_rows < rows; i_rows++) {
    for (std::size_t i_cols = 0; i_cols < cols - 1; i_cols++) {
      if (i_cols == 0) {
        EXPECT_LE(concentrations(i_rows, i_cols), 10);
      } else {
        EXPECT_GE(concentrations(i_rows, i_cols),
                  concentrations(i_rows, i_cols + 1));
      }
    }
  }

  // the values should also be equal from top to bottom
  for (std::size_t i_cols = 0; i_cols < cols; i_cols++) {
    const double &ref = concentrations(0, i_cols);
    for (std::size_t i_rows = 1; i_rows < rows; i_rows++) {
      // check if the values are equal within the epsilon range
      EXPECT_NEAR(ref, concentrations(i_rows, i_cols), epsilon);
    }
  }
}
