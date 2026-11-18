#include <cstddef>
#include <gtest/gtest.h>
#include <tug/tug.hpp>

#define ADVECTION_TEST(x) TEST(Advection, x)

ADVECTION_TEST(LeftToRight) {
  constexpr std::size_t rows = 21;
  constexpr std::size_t cols = 21;

  constexpr double K = 1E-2;
  constexpr double timestep = 500.0;
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
      EXPECT_NEAR(ref, concentrations(i_rows, i_cols), epsilon);
    }
  }
}

ADVECTION_TEST(ExactFrontPropagation1D) {
  constexpr std::size_t cols = 10;
  constexpr double domainL = 10.0; // dx = 1.0
  constexpr double qx = 0.2;       // Darcy velocity = 0.2 m/s
  constexpr double porosity = 0.5; // pore velocity v = qx / porosity = 0.4 m/s
  constexpr double cflDt = 2.5;    // dx / v = 1.0 / 0.4 = 2.5 s
  constexpr int steps = 4;

  tug::RowMajMat<double> hydHeads = tug::RowMajMat<double>::Zero(1, cols);
  tug::RowMajMat<double> concentrations = tug::RowMajMat<double>::Zero(1, cols);

  tug::Velocities<double> velocities(hydHeads);
  velocities.setDomain(domainL);
  velocities.setVelocitiesX(
      tug::RowMajMat<double>::Constant(1, cols + 1, qx));

  tug::Advection advection(concentrations, velocities);
  advection.setPorosity(porosity);
  advection.setTimestep(cflDt);
  advection.setIterations(steps);

  tug::Boundary<double> &bcC = advection.getBoundaryConditions();
  bcC.setBoundarySideConstant(tug::BC_SIDE_LEFT, 1.0);
  bcC.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 0.0);

  advection.run();

  // At Courant number = 1, the step front advances by exactly 1 cell per step
  for (std::size_t j = 0; j < cols; j++) {
    if (j < static_cast<std::size_t>(steps)) {
      EXPECT_NEAR(concentrations(0, j), 1.0, 1E-12);
    } else {
      EXPECT_NEAR(concentrations(0, j), 0.0, 1E-12);
    }
  }
}

ADVECTION_TEST(FourDirections2D) {
  constexpr std::size_t n = 6;
  constexpr double domain = 6.0; // dx = dy = 1.0
  constexpr double q = 0.5;
  constexpr double porosity = 0.5; // v = 1.0 m/s -> cflDt = 1.0 s
  constexpr double dt = 1.0;
  constexpr int steps = 2;

  // 1. Right-to-Left (vx = -q)
  {
    tug::RowMajMat<double> heads = tug::RowMajMat<double>::Zero(n, n);
    tug::RowMajMat<double> conc = tug::RowMajMat<double>::Zero(n, n);
    tug::Velocities<double> velo(heads);
    velo.setDomain(domain, domain);
    velo.setVelocitiesX(tug::RowMajMat<double>::Constant(n, n + 1, -q));
    velo.setVelocitiesY(tug::RowMajMat<double>::Zero(n + 1, n));

    tug::Advection adv(conc, velo);
    adv.setPorosity(porosity);
    adv.setTimestep(dt);
    adv.setIterations(steps);
    adv.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_RIGHT,
                                                        1.0);
    adv.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_LEFT, 0.0);
    adv.run();

    for (std::size_t r = 0; r < n; r++) {
      EXPECT_NEAR(conc(r, n - 1), 1.0, 1E-12);
      EXPECT_NEAR(conc(r, n - 2), 1.0, 1E-12);
      EXPECT_NEAR(conc(r, n - 3), 0.0, 1E-12);
    }
  }

  // 2. Top-to-Bottom (vy = +q)
  {
    tug::RowMajMat<double> heads = tug::RowMajMat<double>::Zero(n, n);
    tug::RowMajMat<double> conc = tug::RowMajMat<double>::Zero(n, n);
    tug::Velocities<double> velo(heads);
    velo.setDomain(domain, domain);
    velo.setVelocitiesX(tug::RowMajMat<double>::Zero(n, n + 1));
    velo.setVelocitiesY(tug::RowMajMat<double>::Constant(n + 1, n, q));

    tug::Advection adv(conc, velo);
    adv.setPorosity(porosity);
    adv.setTimestep(dt);
    adv.setIterations(steps);
    adv.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_TOP, 1.0);
    adv.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_BOTTOM,
                                                        0.0);
    adv.run();

    for (std::size_t c = 0; c < n; c++) {
      EXPECT_NEAR(conc(0, c), 1.0, 1E-12);
      EXPECT_NEAR(conc(1, c), 1.0, 1E-12);
      EXPECT_NEAR(conc(2, c), 0.0, 1E-12);
    }
  }

  // 3. Bottom-to-Top (vy = -q)
  {
    tug::RowMajMat<double> heads = tug::RowMajMat<double>::Zero(n, n);
    tug::RowMajMat<double> conc = tug::RowMajMat<double>::Zero(n, n);
    tug::Velocities<double> velo(heads);
    velo.setDomain(domain, domain);
    velo.setVelocitiesX(tug::RowMajMat<double>::Zero(n, n + 1));
    velo.setVelocitiesY(tug::RowMajMat<double>::Constant(n + 1, n, -q));

    tug::Advection adv(conc, velo);
    adv.setPorosity(porosity);
    adv.setTimestep(dt);
    adv.setIterations(steps);
    adv.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_BOTTOM,
                                                        1.0);
    adv.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_TOP, 0.0);
    adv.run();

    for (std::size_t c = 0; c < n; c++) {
      EXPECT_NEAR(conc(n - 1, c), 1.0, 1E-12);
      EXPECT_NEAR(conc(n - 2, c), 1.0, 1E-12);
      EXPECT_NEAR(conc(n - 3, c), 0.0, 1E-12);
    }
  }
}

ADVECTION_TEST(PorosityScaling) {
  constexpr std::size_t cols = 5;
  constexpr double domainL = 5.0; // dx = 1.0
  constexpr double qx = 0.1;
  constexpr double dt = 1.0; // small enough that innerSteps == 1 for both

  tug::RowMajMat<double> heads1 = tug::RowMajMat<double>::Zero(1, cols);
  tug::RowMajMat<double> conc1 = tug::RowMajMat<double>::Zero(1, cols);
  tug::Velocities<double> velo1(heads1);
  velo1.setDomain(domainL);
  velo1.setVelocitiesX(tug::RowMajMat<double>::Constant(1, cols + 1, qx));

  tug::Advection adv1(conc1, velo1);
  adv1.setPorosity(1.0);
  adv1.setTimestep(dt);
  adv1.setIterations(1);
  adv1.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_LEFT, 1.0);
  adv1.run();

  tug::RowMajMat<double> heads2 = tug::RowMajMat<double>::Zero(1, cols);
  tug::RowMajMat<double> conc2 = tug::RowMajMat<double>::Zero(1, cols);
  tug::Velocities<double> velo2(heads2);
  velo2.setDomain(domainL);
  velo2.setVelocitiesX(tug::RowMajMat<double>::Constant(1, cols + 1, qx));

  tug::Advection adv2(conc2, velo2);
  adv2.setPorosity(0.25);
  adv2.setTimestep(dt);
  adv2.setIterations(1);
  adv2.getBoundaryConditions().setBoundarySideConstant(tug::BC_SIDE_LEFT, 1.0);
  adv2.run();

  // With porosity = 0.25, concentration change in cell 0 after 1 step must be
  // 4 times larger than with porosity = 1.0
  EXPECT_NEAR(conc2(0, 0), 4.0 * conc1(0, 0), 1E-12);
}

ADVECTION_TEST(MassConservationClosedBoundaries) {
  constexpr std::size_t rows = 5;
  constexpr std::size_t cols = 5;

  tug::RowMajMat<double> heads = tug::RowMajMat<double>::Zero(rows, cols);
  tug::RowMajMat<double> conc = tug::RowMajMat<double>::Zero(rows, cols);
  conc(2, 2) = 10.0;
  const double initialMass = conc.sum();

  tug::Velocities<double> velo(heads);
  velo.setDomain(5.0, 5.0);
  velo.setVelocitiesX(tug::RowMajMat<double>::Constant(rows, cols + 1, 0.2));
  velo.setVelocitiesY(tug::RowMajMat<double>::Constant(rows + 1, cols, 0.1));

  // All boundaries are closed by default
  tug::Advection adv(conc, velo);
  adv.setPorosity(0.5);
  adv.setTimestep(5.0);
  adv.setIterations(10);
  adv.run();

  EXPECT_NEAR(conc.sum(), initialMass, 1E-10);
}

ADVECTION_TEST(InnerBoundarySource) {
  constexpr std::size_t rows = 5;
  constexpr std::size_t cols = 5;

  tug::RowMajMat<double> heads = tug::RowMajMat<double>::Zero(rows, cols);
  tug::RowMajMat<double> conc = tug::RowMajMat<double>::Zero(rows, cols);

  tug::Velocities<double> velo(heads);
  velo.setDomain(5.0, 5.0);
  velo.setVelocitiesX(tug::RowMajMat<double>::Constant(rows, cols + 1, 0.5));
  velo.setVelocitiesY(tug::RowMajMat<double>::Zero(rows + 1, cols));

  tug::Advection adv(conc, velo);
  adv.setPorosity(0.5);
  adv.setTimestep(1.0);
  adv.setIterations(3);

  auto &bcC = adv.getBoundaryConditions();
  bcC.setInnerBoundary(2, 2, 5.0);
  bcC.setBoundarySideConstant(tug::BC_SIDE_RIGHT, 0.0);

  adv.run();

  // Source cell must remain constant at 5.0
  EXPECT_DOUBLE_EQ(conc(2, 2), 5.0);
  // Upstream cells (left of source) must remain 0.0
  EXPECT_DOUBLE_EQ(conc(2, 0), 0.0);
  EXPECT_DOUBLE_EQ(conc(2, 1), 0.0);
  // Downstream cells (right of source) must have received concentration
  EXPECT_GT(conc(2, 3), 0.0);
  EXPECT_GT(conc(2, 4), 0.0);
}
