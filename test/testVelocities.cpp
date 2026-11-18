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

VELOCITIES_TEST(HeterogeneousLayeredPermeability1DAnd2D) {
  constexpr std::size_t cols = 10;
  constexpr double domainX = 10.0; // dx = 1.0
  constexpr double K1 = 1E-2;
  constexpr double K2 = 1E-3;
  constexpr double hLeft = 10.0;
  constexpr double hRight = 0.0;

  // Analytical steady-state Darcy flux for two 5m layers in series:
  // q = (hLeft - hRight) / (L1/K1 + L2/K2) = 10 / (5/0.01 + 5/0.001) = 1/550
  const double expectedFlux = (hLeft - hRight) / (5.0 / K1 + 5.0 / K2);

  // 1D test
  {
    tug::RowMajMat<double> heads1D =
        tug::RowMajMat<double>::Constant(1, cols, 5.0);
    tug::RowMajMat<double> perm1D(1, cols);
    for (std::size_t j = 0; j < cols; j++) {
      perm1D(0, j) = (j < cols / 2) ? K1 : K2;
    }

    tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                    tug::HYDRAULIC_RESOLVE::IMPLICIT>
        velo1D(heads1D);
    velo1D.setDomain(domainX);
    velo1D.setPermKX(perm1D);
    velo1D.setEpsilon(1E-11);
    velo1D.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_LEFT,
                                                           hLeft);
    velo1D.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_RIGHT,
                                                           hRight);
    velo1D.run();

    const auto &vx = velo1D.getVelocitiesX();
    for (std::size_t j = 0; j <= cols; j++) {
      EXPECT_NEAR(vx(0, j), expectedFlux, 1E-6);
    }
  }

  // 2D test
  {
    constexpr std::size_t rows = 4;
    tug::RowMajMat<double> heads2D =
        tug::RowMajMat<double>::Constant(rows, cols, 5.0);
    tug::RowMajMat<double> perm2D(rows, cols);
    for (std::size_t i = 0; i < rows; i++) {
      for (std::size_t j = 0; j < cols; j++) {
        perm2D(i, j) = (j < cols / 2) ? K1 : K2;
      }
    }

    tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                    tug::HYDRAULIC_RESOLVE::IMPLICIT>
        velo2D(heads2D);
    velo2D.setDomain(4.0, domainX);
    velo2D.setPermKX(perm2D);
    velo2D.setPermKY(perm2D);
    velo2D.setEpsilon(1E-11);
    velo2D.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_LEFT,
                                                           hLeft);
    velo2D.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_RIGHT,
                                                           hRight);
    velo2D.run();

    const auto &vx = velo2D.getVelocitiesX();
    const auto &vy = velo2D.getVelocitiesY();
    for (std::size_t i = 0; i < rows; i++) {
      for (std::size_t j = 0; j <= cols; j++) {
        EXPECT_NEAR(vx(i, j), expectedFlux, 1E-6);
      }
    }
    for (std::size_t i = 0; i <= rows; i++) {
      for (std::size_t j = 0; j < cols; j++) {
        EXPECT_NEAR(vy(i, j), 0.0, 1E-9);
      }
    }
  }
}

VELOCITIES_TEST(DivergenceFreeAroundInnerSource2D) {
  constexpr std::size_t rows = 15;
  constexpr std::size_t cols = 15;
  constexpr std::size_t center_row = rows / 2;
  constexpr std::size_t center_col = cols / 2;
  constexpr double K = 1E-2;

  tug::RowMajMat<double> hydHeads =
      tug::RowMajMat<double>::Constant(rows, cols, 1.0);
  tug::RowMajMat<double> permK =
      tug::RowMajMat<double>::Constant(rows, cols, K);

  tug::Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
                  tug::HYDRAULIC_RESOLVE::IMPLICIT>
      velo(hydHeads);
  velo.setDomain(15.0, 15.0); // dx = dy = 1.0
  velo.setPermKX(permK);
  velo.setPermKY(permK);
  velo.setEpsilon(1E-10);

  auto &bcH = velo.getBoundaryConditions();
  bcH.setBoundarySideConstant(tug::BC_SIDE_LEFT, 1.0);
  bcH.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 1.0);
  bcH.setBoundarySideConstant(tug::BC_SIDE_TOP, 1.0);
  bcH.setBoundarySideConstant(tug::BC_SIDE_BOTTOM, 1.0);
  bcH.setInnerBoundary(center_row, center_col, 10.0);

  velo.run();

  const auto &vx = velo.getVelocitiesX();
  const auto &vy = velo.getVelocitiesY();

  for (std::size_t i = 0; i < rows; i++) {
    for (std::size_t j = 0; j < cols; j++) {
      const double div = (vx(i, j + 1) - vx(i, j)) + (vy(i + 1, j) - vy(i, j));
      if (i == center_row && j == center_col) {
        // Source cell must have positive net outflow
        EXPECT_GT(div, 0.0);
      } else {
        // Every non-source cell must be divergence-free in steady state
        EXPECT_NEAR(div, 0.0, 1E-6);
      }
    }
  }
}

VELOCITIES_TEST(TransientHydraulicHeadEvolution) {
  constexpr std::size_t cols = 10;
  tug::RowMajMat<double> heads = tug::RowMajMat<double>::Zero(1, cols);
  tug::RowMajMat<double> permK =
      tug::RowMajMat<double>::Constant(1, cols, 1E-2);

  tug::Velocities<double, tug::HYDRAULIC_MODE::TRANSIENT,
                  tug::HYDRAULIC_RESOLVE::IMPLICIT>
      velo(heads);
  velo.setDomain(10.0);
  velo.setPermKX(permK);
  velo.setTimestep(10.0);

  auto &bc = velo.getBoundaryConditions();
  bc.setBoundarySideConstant(tug::BC_SIDE_LEFT, 10.0);
  bc.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 0.0);

  velo.run();
  const double headAfterStep1 = heads(0, 0);
  const double vxInteriorStep1 = velo.getVelocitiesX()(0, 1);
  EXPECT_GT(headAfterStep1, 0.0);
  EXPECT_GT(vxInteriorStep1, 0.0);

  velo.run();
  const double headAfterStep2 = heads(0, 0);
  const double vxInteriorStep2 = velo.getVelocitiesX()(0, 1);
  EXPECT_GT(headAfterStep2, headAfterStep1);
  EXPECT_GT(vxInteriorStep2, vxInteriorStep1);
}
