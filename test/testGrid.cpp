#include <Eigen/Core>
#include <doctest/doctest.h>
#include <tug/Grid.hpp>

TEST_CASE("1D Grid, too small length") {
  int l = 2;
  CHECK_THROWS(Grid(l));
}

TEST_CASE("2D Grid, too small side") {
  int r = 2;
  int c = 4;
  CHECK_THROWS(Grid(r, c));

  r = 4;
  c = 2;
  CHECK_THROWS(Grid(r, c));
}

TEST_CASE("1D Grid") {
  int l = 12;
  Grid grid(l);

  SUBCASE("correct construction") {
    CHECK_EQ(grid.getDim(), 1);
    CHECK_EQ(grid.getLength(), l);
    CHECK_EQ(grid.getCol(), l);
    CHECK_EQ(grid.getRow(), 1);

    CHECK_EQ(grid.getConcentrations().rows(), 1);
    CHECK_EQ(grid.getConcentrations().cols(), l);
    CHECK_EQ(grid.getAlpha().rows(), 1);
    CHECK_EQ(grid.getAlpha().cols(), l);
    CHECK_EQ(grid.getDeltaCol(), 1);

    CHECK_THROWS(grid.getAlphaX());
    CHECK_THROWS(grid.getAlphaY());
    CHECK_THROWS(grid.getDeltaRow());
  }

  SUBCASE("setting concentrations") {
    // correct concentrations matrix
    MatrixXd concentrations = MatrixXd::Constant(1, l, 3);
    CHECK_NOTHROW(grid.setConcentrations(concentrations));

    // false concentrations matrix
    MatrixXd wConcentrations = MatrixXd::Constant(2, l, 4);
    CHECK_THROWS(grid.setConcentrations(wConcentrations));
  }

  SUBCASE("setting alpha") {
    // correct alpha matrix
    MatrixXd alpha = MatrixXd::Constant(1, l, 3);
    CHECK_NOTHROW(grid.setAlpha(alpha));

    CHECK_THROWS(grid.setAlpha(alpha, alpha));

    grid.setAlpha(alpha);
    CHECK_EQ(grid.getAlpha(), alpha);
    CHECK_THROWS(grid.getAlphaX());
    CHECK_THROWS(grid.getAlphaY());

    // false alpha matrix
    MatrixXd wAlpha = MatrixXd::Constant(3, l, 2);
    CHECK_THROWS(grid.setAlpha(wAlpha));
  }

  SUBCASE("setting domain") {
    int d = 8;
    // set 1D domain
    CHECK_NOTHROW(grid.setDomain(d));

    // set 2D domain
    CHECK_THROWS(grid.setDomain(d, d));

    grid.setDomain(d);
    CHECK_EQ(grid.getDeltaCol(), double(d) / double(l));
    CHECK_THROWS(grid.getDeltaRow());

    // set too small domain
    d = 0;
    CHECK_THROWS(grid.setDomain(d));
    d = -2;
    CHECK_THROWS(grid.setDomain(d));
  }
}

TEST_CASE("2D Grid quadratic") {
  int rc = 12;
  Grid grid(rc, rc);

  SUBCASE("correct construction") {
    CHECK_EQ(grid.getDim(), 2);
    CHECK_THROWS(grid.getLength());
    CHECK_EQ(grid.getCol(), rc);
    CHECK_EQ(grid.getRow(), rc);

    CHECK_EQ(grid.getConcentrations().rows(), rc);
    CHECK_EQ(grid.getConcentrations().cols(), rc);
    CHECK_THROWS(grid.getAlpha());

    CHECK_EQ(grid.getAlphaX().rows(), rc);
    CHECK_EQ(grid.getAlphaX().cols(), rc);
    CHECK_EQ(grid.getAlphaY().rows(), rc);
    CHECK_EQ(grid.getAlphaY().cols(), rc);
    CHECK_EQ(grid.getDeltaRow(), 1);
    CHECK_EQ(grid.getDeltaCol(), 1);
  }

  SUBCASE("setting concentrations") {
    // correct concentrations matrix
    MatrixXd concentrations = MatrixXd::Constant(rc, rc, 2);
    CHECK_NOTHROW(grid.setConcentrations(concentrations));

    // false concentrations matrix
    MatrixXd wConcentrations = MatrixXd::Constant(rc, rc + 3, 1);
    CHECK_THROWS(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(rc + 3, rc, 4);
    CHECK_THROWS(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(rc + 2, rc + 4, 6);
    CHECK_THROWS(grid.setConcentrations(wConcentrations));
  }

  SUBCASE("setting alphas") {
    // correct alpha matrices
    MatrixXd alphax = MatrixXd::Constant(rc, rc, 2);
    MatrixXd alphay = MatrixXd::Constant(rc, rc, 4);
    CHECK_NOTHROW(grid.setAlpha(alphax, alphay));

    CHECK_THROWS(grid.setAlpha(alphax));

    grid.setAlpha(alphax, alphay);
    CHECK_EQ(grid.getAlphaX(), alphax);
    CHECK_EQ(grid.getAlphaY(), alphay);
    CHECK_THROWS(grid.getAlpha());

    // false alpha matrices
    alphax = MatrixXd::Constant(rc + 3, rc + 1, 3);
    CHECK_THROWS(grid.setAlpha(alphax, alphay));
    alphay = MatrixXd::Constant(rc + 2, rc + 1, 3);
    CHECK_THROWS(grid.setAlpha(alphax, alphay));
  }

  SUBCASE("setting domain") {
    int dr = 8;
    int dc = 9;

    // set 1D domain
    CHECK_THROWS(grid.setDomain(dr));

    // set 2D domain
    CHECK_NOTHROW(grid.setDomain(dr, dc));

    grid.setDomain(dr, dc);
    CHECK_EQ(grid.getDeltaCol(), double(dc) / double(rc));
    CHECK_EQ(grid.getDeltaRow(), double(dr) / double(rc));

    // set too small domain
    dr = 0;
    CHECK_THROWS(grid.setDomain(dr, dc));
    dr = 8;
    dc = 0;
    CHECK_THROWS(grid.setDomain(dr, dc));
    dr = -2;
    CHECK_THROWS(grid.setDomain(dr, dc));
  }
}

TEST_CASE("2D Grid non-quadratic") {
  int r = 12;
  int c = 15;
  Grid grid(r, c);

  SUBCASE("correct construction") {
    CHECK_EQ(grid.getDim(), 2);
    CHECK_THROWS(grid.getLength());
    CHECK_EQ(grid.getCol(), c);
    CHECK_EQ(grid.getRow(), r);

    CHECK_EQ(grid.getConcentrations().rows(), r);
    CHECK_EQ(grid.getConcentrations().cols(), c);
    CHECK_THROWS(grid.getAlpha());

    CHECK_EQ(grid.getAlphaX().rows(), r);
    CHECK_EQ(grid.getAlphaX().cols(), c);
    CHECK_EQ(grid.getAlphaY().rows(), r);
    CHECK_EQ(grid.getAlphaY().cols(), c);
    CHECK_EQ(grid.getDeltaRow(), 1);
    CHECK_EQ(grid.getDeltaCol(), 1);
  }

  SUBCASE("setting concentrations") {
    // correct concentrations matrix
    MatrixXd concentrations = MatrixXd::Constant(r, c, 2);
    CHECK_NOTHROW(grid.setConcentrations(concentrations));

    // false concentrations matrix
    MatrixXd wConcentrations = MatrixXd::Constant(r, c + 3, 6);
    CHECK_THROWS(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(r + 3, c, 3);
    CHECK_THROWS(grid.setConcentrations(wConcentrations));
    wConcentrations = MatrixXd::Constant(r + 2, c + 4, 2);
    CHECK_THROWS(grid.setConcentrations(wConcentrations));
  }

  SUBCASE("setting alphas") {
    // correct alpha matrices
    MatrixXd alphax = MatrixXd::Constant(r, c, 2);
    MatrixXd alphay = MatrixXd::Constant(r, c, 4);
    CHECK_NOTHROW(grid.setAlpha(alphax, alphay));

    CHECK_THROWS(grid.setAlpha(alphax));

    grid.setAlpha(alphax, alphay);
    CHECK_EQ(grid.getAlphaX(), alphax);
    CHECK_EQ(grid.getAlphaY(), alphay);
    CHECK_THROWS(grid.getAlpha());

    // false alpha matrices
    alphax = MatrixXd::Constant(r + 3, c + 1, 3);
    CHECK_THROWS(grid.setAlpha(alphax, alphay));
    alphay = MatrixXd::Constant(r + 2, c + 1, 5);
    CHECK_THROWS(grid.setAlpha(alphax, alphay));
  }

  SUBCASE("setting domain") {
    int dr = 8;
    int dc = 9;

    // set 1D domain
    CHECK_THROWS(grid.setDomain(dr));

    // set 2D domain
    CHECK_NOTHROW(grid.setDomain(dr, dc));

    grid.setDomain(dr, dc);
    CHECK_EQ(grid.getDeltaCol(), double(dc) / double(c));
    CHECK_EQ(grid.getDeltaRow(), double(dr) / double(r));

    // set too small domain
    dr = 0;
    CHECK_THROWS(grid.setDomain(dr, dc));
    dr = 8;
    dc = -1;
    CHECK_THROWS(grid.setDomain(dr, dc));
    dr = -2;
    CHECK_THROWS(grid.setDomain(dr, dc));
  }
}