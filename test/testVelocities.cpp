#include <cstddef>
#include <gtest/gtest.h>
#include <tug/Advection/Velocities.hpp>
#include <tug/Boundary.hpp>
#include <tug/Core/Matrix.hpp>

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

  tug::RowMajMat<double> permK =
      tug::RowMajMat<double>::Constant(rows, cols, K);

  tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                  tug::HYDRAULIC_RESOLVE::EXPLICIT>
      velo(hydHeads);

  velo.setDomain(100, 100);
  velo.setPermKX(permK);
  velo.setPermKY(permK);

  tug::Boundary<double> &bcH = velo.getBoundaryConditions();

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

VELOCITIES_TEST(UniformFlow1D) {
  constexpr std::size_t cols = 20;
  constexpr double domainL = 100.0;
  constexpr double K = 1E-2;
  constexpr double hLeft = 10.0;
  constexpr double hRight = 0.0;

  tug::RowMajMat<double> hydHeads = tug::RowMajMat<double>::Constant(1, cols, 5.0);
  tug::RowMajMat<double> permK = tug::RowMajMat<double>::Constant(1, cols, K);

  tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                  tug::HYDRAULIC_RESOLVE::IMPLICIT>
      velo(hydHeads);

  velo.setDomain(domainL);
  velo.setPermKX(permK);
  velo.setEpsilon(1E-10);

  tug::Boundary<double> &bc = velo.getBoundaryConditions();
  bc.setBoundarySideConstant(tug::BC_SIDE_LEFT, hLeft);
  bc.setBoundarySideConstant(tug::BC_SIDE_RIGHT, hRight);

  velo.run();

  const double expectedVelocity = K * (hLeft - hRight) / domainL; // 1E-3
  const auto &vx = velo.getVelocitiesX();

  ASSERT_EQ(vx.rows(), 1);
  ASSERT_EQ(vx.cols(), cols + 1);

  for (std::size_t j = 0; j <= cols; j++) {
    EXPECT_NEAR(vx(0, j), expectedVelocity, 1E-5);
  }
}

VELOCITIES_TEST(UniformFlow2DExplicitAndImplicit) {
  constexpr std::size_t rows = 5;
  constexpr std::size_t cols = 10;
  constexpr double domainY = 50.0;
  constexpr double domainX = 100.0;
  constexpr double K = 1E-2;
  constexpr double hLeft = 10.0;
  constexpr double hRight = 0.0;

  tug::RowMajMat<double> hydHeadsExp =
      tug::RowMajMat<double>::Constant(rows, cols, 5.0);
  tug::RowMajMat<double> hydHeadsImp =
      tug::RowMajMat<double>::Constant(rows, cols, 5.0);
  tug::RowMajMat<double> permK =
      tug::RowMajMat<double>::Constant(rows, cols, K);

  tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                  tug::HYDRAULIC_RESOLVE::EXPLICIT>
      veloExp(hydHeadsExp);
  veloExp.setDomain(domainY, domainX);
  veloExp.setPermKX(permK);
  veloExp.setPermKY(permK);
  veloExp.setEpsilon(1E-9);
  veloExp.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_LEFT,
                                                          hLeft);
  veloExp.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_RIGHT,
                                                          hRight);
  veloExp.run();

  tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                  tug::HYDRAULIC_RESOLVE::IMPLICIT>
      veloImp(hydHeadsImp);
  veloImp.setDomain(domainY, domainX);
  veloImp.setPermKX(permK);
  veloImp.setPermKY(permK);
  veloImp.setEpsilon(1E-9);
  veloImp.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_LEFT,
                                                          hLeft);
  veloImp.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_RIGHT,
                                                          hRight);
  veloImp.run();

  const double expectedVx = K * (hLeft - hRight) / domainX; // 1E-3

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = 0; j <= cols; j++) {
      EXPECT_NEAR(veloExp.getVelocitiesX()(i, j), expectedVx, 1E-5);
      EXPECT_NEAR(veloImp.getVelocitiesX()(i, j), expectedVx, 1E-5);
    }
  }

  for (std::size_t i = 0; i <= rows; i++) {
    for (std::size_t j = 0; j < cols; j++) {
      EXPECT_NEAR(veloExp.getVelocitiesY()(i, j), 0.0, 1E-8);
      EXPECT_NEAR(veloImp.getVelocitiesY()(i, j), 0.0, 1E-8);
    }
  }
}

VELOCITIES_TEST(CustomVelocitiesSetter) {
  constexpr std::size_t rows = 4;
  constexpr std::size_t cols = 6;

  tug::RowMajMat<double> hydHeads =
      tug::RowMajMat<double>::Zero(rows, cols);
  tug::Velocities<double> velo(hydHeads);

  tug::RowMajMat<double> customVx =
      tug::RowMajMat<double>::Constant(rows, cols + 1, 0.25);
  tug::RowMajMat<double> customVy =
      tug::RowMajMat<double>::Constant(rows + 1, cols, -0.5);

  velo.setVelocitiesX(customVx);
  velo.setVelocitiesY(customVy);
  velo.run();

  EXPECT_DOUBLE_EQ(velo.getVelocitiesX().minCoeff(), 0.25);
  EXPECT_DOUBLE_EQ(velo.getVelocitiesX().maxCoeff(), 0.25);
  EXPECT_DOUBLE_EQ(velo.getVelocitiesY().minCoeff(), -0.5);
  EXPECT_DOUBLE_EQ(velo.getVelocitiesY().maxCoeff(), -0.5);
}
