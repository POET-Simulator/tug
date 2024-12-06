#include <Eigen/Core>
#include <tug/Grid.hpp>
#include <vector>

#include <gtest/gtest.h>

using namespace Eigen;
using namespace std;
using namespace tug;

#define GRID_TEST(x) TEST(Grid, x)

GRID_TEST(InvalidConstructorParams) {
  EXPECT_ANY_THROW(Grid64(2));
  EXPECT_ANY_THROW(Grid64(1, 4));
  EXPECT_ANY_THROW(Grid64(4, 1));
}

GRID_TEST(Grid64OneDimensional) {
  int l = 12;
  Grid64 grid(l);

  {
    EXPECT_EQ(grid.getDim(), 1);
    EXPECT_EQ(grid.getLength(), l);
    EXPECT_EQ(grid.getCol(), l);
    EXPECT_EQ(grid.getRow(), 1);

    EXPECT_EQ(grid.getConcentrations().rows(), 1);
    EXPECT_EQ(grid.getConcentrations().cols(), l);
    EXPECT_EQ(grid.getAlpha().rows(), 1);
    EXPECT_EQ(grid.getAlpha().cols(), l);
    EXPECT_EQ(grid.getDeltaCol(), 1);

    EXPECT_ANY_THROW(grid.getAlphaX());
    EXPECT_ANY_THROW(grid.getAlphaY());
    EXPECT_ANY_THROW(grid.getDeltaRow());
  }

  {
    // correct concentrations matrix
    MatrixXd concentrations = MatrixXd::Constant(1, l, 3);
    EXPECT_NO_THROW(grid.setConcentrations(concentrations));

    // false concentrations matrix
    MatrixXd wConcentrations = MatrixXd::Constant(2, l, 4);
    EXPECT_ANY_THROW(grid.setConcentrations(wConcentrations));
  }

  {
    // correct alpha matrix
    MatrixXd alpha = MatrixXd::Constant(1, l, 3);
    EXPECT_NO_THROW(grid.setAlpha(alpha));

    EXPECT_ANY_THROW(grid.setAlpha(alpha, alpha));

    grid.setAlpha(alpha);
    EXPECT_EQ(grid.getAlpha(), alpha);
    EXPECT_ANY_THROW(grid.getAlphaX());
    EXPECT_ANY_THROW(grid.getAlphaY());

    // false alpha matrix
    MatrixXd wAlpha = MatrixXd::Constant(3, l, 2);
    EXPECT_ANY_THROW(grid.setAlpha(wAlpha));
  }

  {
    int d = 8;
    // set 1D domain
    EXPECT_NO_THROW(grid.setDomain(d));

    // set 2D domain
    EXPECT_ANY_THROW(grid.setDomain(d, d));

    grid.setDomain(d);
    EXPECT_DOUBLE_EQ(grid.getDeltaCol(), double(d) / double(l));
    EXPECT_ANY_THROW(grid.getDeltaRow());

    // set too small domain
    d = 0;
    EXPECT_ANY_THROW(grid.setDomain(d));
    d = -2;
    EXPECT_ANY_THROW(grid.setDomain(d));
  }
}

GRID_TEST(Grid64Quadratic) {
  int rc = 12;
  Grid64 grid(rc, rc);

  {
    EXPECT_EQ(grid.getDim(), 2);
    EXPECT_ANY_THROW(grid.getLength());
    EXPECT_EQ(grid.getCol(), rc);
    EXPECT_EQ(grid.getRow(), rc);

    EXPECT_EQ(grid.getConcentrations().rows(), rc);
    EXPECT_EQ(grid.getConcentrations().cols(), rc);
    EXPECT_ANY_THROW(grid.getAlpha());

    EXPECT_EQ(grid.getAlphaX().rows(), rc);
    EXPECT_EQ(grid.getAlphaX().cols(), rc);
    EXPECT_EQ(grid.getAlphaY().rows(), rc);
    EXPECT_EQ(grid.getAlphaY().cols(), rc);
    EXPECT_EQ(grid.getDeltaRow(), 1);
    EXPECT_EQ(grid.getDeltaCol(), 1);
  }

  {
    // correct concentrations matrix
    MatrixXd concentrations = MatrixXd::Constant(rc, rc, 2);
    EXPECT_NO_THROW(grid.setConcentrations(concentrations));

    // false concentrations matrix
    MatrixXd wConcentrations = MatrixXd::Constant(rc, rc + 3, 1);
    EXPECT_ANY_THROW(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(rc + 3, rc, 4);
    EXPECT_ANY_THROW(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(rc + 2, rc + 4, 6);
    EXPECT_ANY_THROW(grid.setConcentrations(wConcentrations));
  }

  {
    // correct alpha matrices
    MatrixXd alphax = MatrixXd::Constant(rc, rc, 2);
    MatrixXd alphay = MatrixXd::Constant(rc, rc, 4);
    EXPECT_NO_THROW(grid.setAlpha(alphax, alphay));

    EXPECT_ANY_THROW(grid.setAlpha(alphax));

    grid.setAlpha(alphax, alphay);
    EXPECT_EQ(grid.getAlphaX(), alphax);
    EXPECT_EQ(grid.getAlphaY(), alphay);
    EXPECT_ANY_THROW(grid.getAlpha());

    // false alpha matrices
    alphax = MatrixXd::Constant(rc + 3, rc + 1, 3);
    EXPECT_ANY_THROW(grid.setAlpha(alphax, alphay));
    alphay = MatrixXd::Constant(rc + 2, rc + 1, 3);
    EXPECT_ANY_THROW(grid.setAlpha(alphax, alphay));
  }

  {
    int dr = 8;
    int dc = 9;

    // set 1D domain
    EXPECT_ANY_THROW(grid.setDomain(dr));

    // set 2D domain
    EXPECT_NO_THROW(grid.setDomain(dr, dc));

    grid.setDomain(dr, dc);
    EXPECT_DOUBLE_EQ(grid.getDeltaCol(), double(dc) / double(rc));
    EXPECT_DOUBLE_EQ(grid.getDeltaRow(), double(dr) / double(rc));

    // set too small domain
    dr = 0;
    EXPECT_ANY_THROW(grid.setDomain(dr, dc));
    dr = 8;
    dc = 0;
    EXPECT_ANY_THROW(grid.setDomain(dr, dc));
    dr = -2;
    EXPECT_ANY_THROW(grid.setDomain(dr, dc));
  }
}

GRID_TEST(Grid64NonQuadratic) {
  int r = 12;
  int c = 15;
  Grid64 grid(r, c);

  {
    EXPECT_EQ(grid.getDim(), 2);
    EXPECT_ANY_THROW(grid.getLength());
    EXPECT_EQ(grid.getCol(), c);
    EXPECT_EQ(grid.getRow(), r);

    EXPECT_EQ(grid.getConcentrations().rows(), r);
    EXPECT_EQ(grid.getConcentrations().cols(), c);
    EXPECT_ANY_THROW(grid.getAlpha());

    EXPECT_EQ(grid.getAlphaX().rows(), r);
    EXPECT_EQ(grid.getAlphaX().cols(), c);
    EXPECT_EQ(grid.getAlphaY().rows(), r);
    EXPECT_EQ(grid.getAlphaY().cols(), c);
    EXPECT_EQ(grid.getDeltaRow(), 1);
    EXPECT_EQ(grid.getDeltaCol(), 1);
  }

  {
    // correct concentrations matrix
    MatrixXd concentrations = MatrixXd::Constant(r, c, 2);
    EXPECT_NO_THROW(grid.setConcentrations(concentrations));

    // false concentrations matrix
    MatrixXd wConcentrations = MatrixXd::Constant(r, c + 3, 6);
    EXPECT_ANY_THROW(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(r + 3, c, 3);
    EXPECT_ANY_THROW(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(r + 2, c + 4, 2);
    EXPECT_ANY_THROW(grid.setConcentrations(wConcentrations));
  }

  {
    // correct alpha matrices
    MatrixXd alphax = MatrixXd::Constant(r, c, 2);
    MatrixXd alphay = MatrixXd::Constant(r, c, 4);
    EXPECT_NO_THROW(grid.setAlpha(alphax, alphay));

    EXPECT_ANY_THROW(grid.setAlpha(alphax));

    grid.setAlpha(alphax, alphay);
    EXPECT_EQ(grid.getAlphaX(), alphax);
    EXPECT_EQ(grid.getAlphaY(), alphay);
    EXPECT_ANY_THROW(grid.getAlpha());

    // false alpha matrices
    alphax = MatrixXd::Constant(r + 3, c + 1, 3);
    EXPECT_ANY_THROW(grid.setAlpha(alphax, alphay));
    alphay = MatrixXd::Constant(r + 2, c + 1, 5);
    EXPECT_ANY_THROW(grid.setAlpha(alphax, alphay));
  }

  {
    int dr = 8;
    int dc = 9;

    // set 1D domain
    EXPECT_ANY_THROW(grid.setDomain(dr));

    // set 2D domain
    EXPECT_NO_THROW(grid.setDomain(dr, dc));

    grid.setDomain(dr, dc);
    EXPECT_DOUBLE_EQ(grid.getDeltaCol(), double(dc) / double(c));
    EXPECT_DOUBLE_EQ(grid.getDeltaRow(), double(dr) / double(r));

    // set too small domain
    dr = 0;
    EXPECT_ANY_THROW(grid.setDomain(dr, dc));
    dr = 8;
    dc = -1;
    EXPECT_ANY_THROW(grid.setDomain(dr, dc));
    dr = -2;
    EXPECT_ANY_THROW(grid.setDomain(dr, dc));
  }

  {
    std::vector<double> concentrations(r * c);

    for (int i = 0; i < r * c; i++) {
      concentrations[i] = i;
    }

    grid.setConcentrations(concentrations.data());

    EXPECT_DOUBLE_EQ(grid.getConcentrations()(0, 0), 0);
    EXPECT_DOUBLE_EQ(grid.getConcentrations()(0, 1), 1);
    EXPECT_DOUBLE_EQ(grid.getConcentrations()(1, 0), c);
  }
}
