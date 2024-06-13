#include "gtest/gtest.h"
#include <Eigen/Core>
#include <tug/Grid.hpp>

#include <gtest/gtest.h>

using namespace Eigen;
using namespace std;
using namespace tug;

#define GRID_TEST(x) TEST(Grid, x)

GRID_TEST(Grid64OneDimensional) {
  int l = 12;
  Eigen::VectorXd conc(l);
  Grid64 grid(conc);
  grid.initAlpha();

  {
    EXPECT_EQ(grid.getDim(), 1);
    EXPECT_EQ(grid.getCol(), l);
    EXPECT_EQ(grid.getCol(), l);
    EXPECT_EQ(grid.getRow(), 1);

    EXPECT_EQ(grid.getConcentrations().rows(), 1);
    EXPECT_EQ(grid.getConcentrations().cols(), l);
    EXPECT_EQ(grid.getAlphaX().rows(), 1);
    EXPECT_EQ(grid.getAlphaX().cols(), l);
    EXPECT_EQ(grid.getDeltaCol(), 1);

    EXPECT_DEATH(grid.getAlphaY(), ".* no alphaY!.*");
    EXPECT_DEATH(grid.getDeltaRow(), ".* not two dimensional, .*");
  }

  {
    // correct alpha matrix
    MatrixXd alpha = MatrixXd::Constant(1, l, 3);
    EXPECT_NO_THROW(grid.setAlpha(alpha));

    EXPECT_DEATH(grid.setAlpha(alpha, alpha), ".* is not two dimensional, .*");

    grid.setAlpha(alpha);
    EXPECT_EQ(grid.getAlphaX(), alpha);
    EXPECT_NO_THROW(grid.getAlphaX());
    EXPECT_DEATH(grid.getAlphaY(), ".* no alphaY!.*");

    // false alpha matrix
    MatrixXd wAlpha = MatrixXd::Constant(3, l, 2);
    EXPECT_DEATH(grid.setAlpha(wAlpha), ".* mismatch with Grid dimensions!.*");
  }

  {
    int d = 8;
    // set 1D domain
    EXPECT_NO_THROW(grid.setDomain(d));

    // set 2D domain
    EXPECT_DEATH(grid.setDomain(d, d), ".* not two dimensional, .*");

    grid.setDomain(d);
    EXPECT_DOUBLE_EQ(grid.getDeltaCol(), double(d) / double(l));
    EXPECT_DEATH(grid.getDeltaRow(), ".* not two dimensional, .*");

    // set too small domain
    EXPECT_DEATH(grid.setDomain(-2), "Given domain length .*");
  }
}

GRID_TEST(Grid64Quadratic) {
  int rc = 12;
  Eigen::MatrixXd conc(rc, rc);
  Grid64 grid(conc);
  grid.initAlpha();

  {
    EXPECT_EQ(grid.getDim(), 2);
    EXPECT_EQ(grid.getCol(), rc);
    EXPECT_EQ(grid.getRow(), rc);

    EXPECT_EQ(grid.getConcentrations().rows(), rc);
    EXPECT_EQ(grid.getConcentrations().cols(), rc);

    EXPECT_EQ(grid.getAlphaX().rows(), rc);
    EXPECT_EQ(grid.getAlphaX().cols(), rc);
    EXPECT_EQ(grid.getAlphaY().rows(), rc);
    EXPECT_EQ(grid.getAlphaY().cols(), rc);
    EXPECT_EQ(grid.getDeltaRow(), 1);
    EXPECT_EQ(grid.getDeltaCol(), 1);
  }
  {
    // correct alpha matrices
    MatrixXd alphax = MatrixXd::Constant(rc, rc, 2);
    MatrixXd alphay = MatrixXd::Constant(rc, rc, 4);
    EXPECT_NO_THROW(grid.setAlpha(alphax, alphay));

    EXPECT_DEATH(grid.setAlpha(alphax), ".* 2D setter function!.*");

    grid.setAlpha(alphax, alphay);
    EXPECT_EQ(grid.getAlphaX(), alphax);
    EXPECT_EQ(grid.getAlphaY(), alphay);

    // false alpha matrices
    alphax = MatrixXd::Constant(rc + 3, rc + 1, 3);
    EXPECT_DEATH(grid.setAlpha(alphax, alphay),
                 ".*has wrong number of rows!.*");
    alphay = MatrixXd::Constant(rc + 2, rc + 1, 3);
    EXPECT_DEATH(grid.setAlpha(alphax, alphay),
                 ".*has wrong number of rows!.*");
  }

  {
    int dr = 8;
    int dc = 9;

    // set 1D domain
    EXPECT_DEATH(grid.setDomain(dr), ".* 2D domain setter!.*");

    // set 2D domain
    EXPECT_NO_THROW(grid.setDomain(dr, dc));

    grid.setDomain(dr, dc);
    EXPECT_DOUBLE_EQ(grid.getDeltaCol(), double(dc) / double(rc));
    EXPECT_DOUBLE_EQ(grid.getDeltaRow(), double(dr) / double(rc));

    // set too small domain
    dr = 0;
    EXPECT_DEATH(grid.setDomain(dr, dc), ".* not positive!.*");
    dr = 8;
    dc = 0;
    EXPECT_DEATH(grid.setDomain(dr, dc), ".* not positive!.*");
    dr = -2;
    EXPECT_DEATH(grid.setDomain(dr, dc), ".* not positive!.*");
  }
}

GRID_TEST(Grid64NonQuadratic) {
  int r = 12;
  int c = 15;
  Eigen::MatrixXd conc(r, c);
  Grid64 grid(conc);
  grid.initAlpha();

  {
    EXPECT_EQ(grid.getDim(), 2);
    EXPECT_EQ(grid.getCol(), c);
    EXPECT_EQ(grid.getRow(), r);

    EXPECT_EQ(grid.getConcentrations().rows(), r);
    EXPECT_EQ(grid.getConcentrations().cols(), c);

    EXPECT_EQ(grid.getAlphaX().rows(), r);
    EXPECT_EQ(grid.getAlphaX().cols(), c);
    EXPECT_EQ(grid.getAlphaY().rows(), r);
    EXPECT_EQ(grid.getAlphaY().cols(), c);
    EXPECT_EQ(grid.getDeltaRow(), 1);
    EXPECT_EQ(grid.getDeltaCol(), 1);
  }

  {
    // correct alpha matrices
    MatrixXd alphax = MatrixXd::Constant(r, c, 2);
    MatrixXd alphay = MatrixXd::Constant(r, c, 4);
    EXPECT_NO_THROW(grid.setAlpha(alphax, alphay));

    grid.setAlpha(alphax, alphay);
    EXPECT_EQ(grid.getAlphaX(), alphax);
    EXPECT_EQ(grid.getAlphaY(), alphay);
  }

  {
    int dr = 8;
    int dc = 9;

    // set 2D domain
    EXPECT_NO_THROW(grid.setDomain(dr, dc));

    grid.setDomain(dr, dc);
    EXPECT_EQ(grid.getDeltaCol(), double(dc) / double(c));
    EXPECT_EQ(grid.getDeltaRow(), double(dr) / double(r));
  }

  {
    int r = 4;
    int c = 5;
    std::vector<double> concentrations(r * c);

    for (int i = 0; i < r * c; i++) {
      concentrations[i] = i;
    }
    Grid64 grid(concentrations.data(), r, c);
    grid.initAlpha();

    {
      EXPECT_EQ(grid.getDim(), 2);
      EXPECT_EQ(grid.getCol(), c);
      EXPECT_EQ(grid.getRow(), r);

      EXPECT_EQ(grid.getConcentrations().rows(), r);
      EXPECT_EQ(grid.getConcentrations().cols(), c);

      EXPECT_EQ(grid.getAlphaX().rows(), r);
      EXPECT_EQ(grid.getAlphaX().cols(), c);
      EXPECT_EQ(grid.getAlphaY().rows(), r);
      EXPECT_EQ(grid.getAlphaY().cols(), c);
      EXPECT_EQ(grid.getDeltaRow(), 1);
      EXPECT_EQ(grid.getDeltaCol(), 1);
    }

    {
      // correct alpha matrices
      MatrixXd alphax = MatrixXd::Constant(r, c, 2);
      MatrixXd alphay = MatrixXd::Constant(r, c, 4);
      EXPECT_NO_THROW(grid.setAlpha(alphax, alphay));

      EXPECT_DEATH(grid.setAlpha(alphax), ".* 2D setter function!.*");

      grid.setAlpha(alphax, alphay);
      EXPECT_EQ(grid.getAlphaX(), alphax);
      EXPECT_EQ(grid.getAlphaY(), alphay);

      // false alpha matrices
      alphax = MatrixXd::Constant(r + 3, c + 1, 3);
      EXPECT_DEATH(grid.setAlpha(alphax, alphay),
                   ".*has wrong number of rows!.*");
      alphay = MatrixXd::Constant(r + 2, c + 1, 5);
      EXPECT_DEATH(grid.setAlpha(alphax, alphay),
                   ".*has wrong number of rows!.*");

      {
        int dr = 8;
        int dc = 9;

        // set 1D domain
        EXPECT_DEATH(grid.setDomain(dr), ".* 2D domain setter!.*");

        // set 2D domain
        EXPECT_NO_THROW(grid.setDomain(dr, dc));

        grid.setDomain(dr, dc);
        EXPECT_DOUBLE_EQ(grid.getDeltaCol(), double(dc) / double(c));
        EXPECT_DOUBLE_EQ(grid.getDeltaRow(), double(dr) / double(r));
      }

      {
        auto &concentrations = grid.getConcentrations();

        for (int i = 0; i < r; i++) {
          for (int j = 0; j < c; j++) {
            concentrations(i, j) = i * c + j;
          }
        }

        EXPECT_DOUBLE_EQ(grid.getConcentrations()(0, 0), 0);
        EXPECT_DOUBLE_EQ(grid.getConcentrations()(0, 1), 1);
        EXPECT_DOUBLE_EQ(grid.getConcentrations()(1, 0), c);

        EXPECT_DOUBLE_EQ(grid.getConcentrations()(0, 0), 0);
        EXPECT_DOUBLE_EQ(grid.getConcentrations()(0, 1), 1);
        EXPECT_DOUBLE_EQ(grid.getConcentrations()(1, 0), c);
        EXPECT_DOUBLE_EQ(grid.getConcentrations()(2, 1), 2 * c + 1);
      }
    }
  }
}