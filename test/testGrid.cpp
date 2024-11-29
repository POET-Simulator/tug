#include "tug/UniformGrid.hpp"
#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>
#include <doctest/doctest.h>
#include <tug/UniformGrid.hpp>
#include <vector>

using namespace Eigen;
using namespace std;
using namespace tug;

TEST_CASE("1D Grid64") {
  int l = 12;
  Eigen::VectorXd conc(l);
  UniformGrid64 grid(l, l, 1);

  SUBCASE("correct construction") {
    CHECK_EQ(grid.getDim(), 1);
    CHECK_EQ(grid.getCol(), l);
    CHECK_EQ(grid.getCol(), l);
    CHECK_EQ(grid.getRow(), 1);
    CHECK_EQ(grid.getDeltaCol(), 1);

    // CHECK_EQ(grid.getConcentrations().rows(), 1);
    // CHECK_EQ(grid.getConcentrations().cols(), l);
    // CHECK_EQ(grid.getAlphaX().rows(), 1);
    // CHECK_EQ(grid.getAlphaX().cols(), l);
  }

  // SUBCASE("setting alpha") {
  //   // correct alpha matrix
  //   MatrixXd alpha = MatrixXd::Constant(1, l, 3);
  //   CHECK_NOTHROW(grid.setAlpha(alpha));

  //   grid.setAlpha(alpha);
  //   CHECK_EQ(grid.getAlphaX(), alpha);
  // }

  SUBCASE("setting domain") {
    int d = 8;
    // set 1D domain
    CHECK_NOTHROW(grid.setDomain(d));

    grid.setDomain(d);
    CHECK_EQ(grid.getDeltaCol(), double(d) / double(l));
  }
}

// TEST_CASE("2D Grid64 quadratic") {
//   int rc = 12;
//   Eigen::MatrixXd conc(rc, rc);
//   Grid64 grid(conc);
//   grid.initAlpha();

//   SUBCASE("correct construction") {
//     CHECK_EQ(grid.getDim(), 2);
//     CHECK_EQ(grid.getCol(), rc);
//     CHECK_EQ(grid.getRow(), rc);

//     CHECK_EQ(grid.getConcentrations().rows(), rc);
//     CHECK_EQ(grid.getConcentrations().cols(), rc);

//     CHECK_EQ(grid.getAlphaX().rows(), rc);
//     CHECK_EQ(grid.getAlphaX().cols(), rc);
//     CHECK_EQ(grid.getAlphaY().rows(), rc);
//     CHECK_EQ(grid.getAlphaY().cols(), rc);
//     CHECK_EQ(grid.getDeltaRow(), 1);
//     CHECK_EQ(grid.getDeltaCol(), 1);
//   }

//   SUBCASE("setting alphas") {
//     // correct alpha matrices
//     MatrixXd alphax = MatrixXd::Constant(rc, rc, 2);
//     MatrixXd alphay = MatrixXd::Constant(rc, rc, 4);
//     CHECK_NOTHROW(grid.setAlpha(alphax, alphay));

//     grid.setAlpha(alphax, alphay);
//     CHECK_EQ(grid.getAlphaX(), alphax);
//     CHECK_EQ(grid.getAlphaY(), alphay);
//   }

//   SUBCASE("setting domain") {
//     int dr = 8;
//     int dc = 9;

//     // set 2D domain
//     CHECK_NOTHROW(grid.setDomain(dr, dc));

//     grid.setDomain(dr, dc);
//     CHECK_EQ(grid.getDeltaCol(), double(dc) / double(rc));
//     CHECK_EQ(grid.getDeltaRow(), double(dr) / double(rc));
//   }
// }

// TEST_CASE("2D Grid64 non-quadratic") {
//   int r = 12;
//   int c = 15;
//   Eigen::MatrixXd conc(r, c);
//   Grid64 grid(conc);
//   grid.initAlpha();

//   SUBCASE("correct construction") {
//     CHECK_EQ(grid.getDim(), 2);
//     CHECK_EQ(grid.getCol(), c);
//     CHECK_EQ(grid.getRow(), r);

//     CHECK_EQ(grid.getConcentrations().rows(), r);
//     CHECK_EQ(grid.getConcentrations().cols(), c);

//     CHECK_EQ(grid.getAlphaX().rows(), r);
//     CHECK_EQ(grid.getAlphaX().cols(), c);
//     CHECK_EQ(grid.getAlphaY().rows(), r);
//     CHECK_EQ(grid.getAlphaY().cols(), c);
//     CHECK_EQ(grid.getDeltaRow(), 1);
//     CHECK_EQ(grid.getDeltaCol(), 1);
//   }

//   SUBCASE("setting alphas") {
//     // correct alpha matrices
//     MatrixXd alphax = MatrixXd::Constant(r, c, 2);
//     MatrixXd alphay = MatrixXd::Constant(r, c, 4);
//     CHECK_NOTHROW(grid.setAlpha(alphax, alphay));

//     grid.setAlpha(alphax, alphay);
//     CHECK_EQ(grid.getAlphaX(), alphax);
//     CHECK_EQ(grid.getAlphaY(), alphay);
//   }

//   SUBCASE("setting domain") {
//     int dr = 8;
//     int dc = 9;

//     // set 2D domain
//     CHECK_NOTHROW(grid.setDomain(dr, dc));

//     grid.setDomain(dr, dc);
//     CHECK_EQ(grid.getDeltaCol(), double(dc) / double(c));
//     CHECK_EQ(grid.getDeltaRow(), double(dr) / double(r));
//   }
// }

// TEST_CASE("2D Grid64 non-quadratic from pointer") {
//   int r = 4;
//   int c = 5;
//   std::vector<double> concentrations(r * c);

//   for (int i = 0; i < r * c; i++) {
//     concentrations[i] = i;
//   }
//   Grid64 grid(concentrations.data(), r, c);
//   grid.initAlpha();

//   SUBCASE("correct construction") {
//     CHECK_EQ(grid.getDim(), 2);
//     CHECK_EQ(grid.getCol(), c);
//     CHECK_EQ(grid.getRow(), r);

//     CHECK_EQ(grid.getConcentrations().rows(), r);
//     CHECK_EQ(grid.getConcentrations().cols(), c);

//     CHECK_EQ(grid.getAlphaX().rows(), r);
//     CHECK_EQ(grid.getAlphaX().cols(), c);
//     CHECK_EQ(grid.getAlphaY().rows(), r);
//     CHECK_EQ(grid.getAlphaY().cols(), c);
//     CHECK_EQ(grid.getDeltaRow(), 1);
//     CHECK_EQ(grid.getDeltaCol(), 1);
//   }

//   SUBCASE("setting alphas") {
//     // correct alpha matrices
//     MatrixXd alphax = MatrixXd::Constant(r, c, 2);
//     MatrixXd alphay = MatrixXd::Constant(r, c, 4);
//     CHECK_NOTHROW(grid.setAlpha(alphax, alphay));

//     grid.setAlpha(alphax, alphay);
//     CHECK_EQ(grid.getAlphaX(), alphax);
//     CHECK_EQ(grid.getAlphaY(), alphay);
//   }

//   SUBCASE("setting domain") {
//     int dr = 8;
//     int dc = 9;

//     // set 2D domain
//     CHECK_NOTHROW(grid.setDomain(dr, dc));

//     grid.setDomain(dr, dc);
//     CHECK_EQ(grid.getDeltaCol(), double(dc) / double(c));
//     CHECK_EQ(grid.getDeltaRow(), double(dr) / double(r));
//   }

//   SUBCASE("correct values") {
//     CHECK_EQ(grid.getConcentrations()(0, 0), 0);
//     CHECK_EQ(grid.getConcentrations()(0, 1), 1);
//     CHECK_EQ(grid.getConcentrations()(1, 0), c);
//     CHECK_EQ(grid.getConcentrations()(2, 1), 2 * c + 1);
//   }
// }