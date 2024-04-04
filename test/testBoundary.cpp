#include <doctest/doctest.h>
#include <iostream>
#include <stdio.h>
#include <string>
#include <tug/Boundary.hpp>
#include <typeinfo>
#include <utility>
#include <vector>

using namespace std;
using namespace tug;

TEST_CASE("BoundaryElement") {

  SUBCASE("Closed case") {
    BoundaryElement boundaryElementClosed = BoundaryElement<double>();
    CHECK_NOTHROW(BoundaryElement<double>());
    CHECK_EQ(boundaryElementClosed.getType(), BC_TYPE_CLOSED);
    CHECK_EQ(boundaryElementClosed.getValue(), -1);
    CHECK_THROWS(boundaryElementClosed.setValue(0.2));
  }

  SUBCASE("Constant case") {
    BoundaryElement boundaryElementConstant = BoundaryElement(0.1);
    CHECK_NOTHROW(BoundaryElement(0.1));
    CHECK_EQ(boundaryElementConstant.getType(), BC_TYPE_CONSTANT);
    CHECK_EQ(boundaryElementConstant.getValue(), 0.1);
    CHECK_NOTHROW(boundaryElementConstant.setValue(0.2));
    CHECK_EQ(boundaryElementConstant.getValue(), 0.2);
  }
}

TEST_CASE("Boundary Class") {
  Grid grid1D = Grid64(10);
  Grid grid2D = Grid64(10, 12);
  Boundary boundary1D = Boundary(grid1D);
  Boundary boundary2D = Boundary(grid2D);
  vector<BoundaryElement<double>> boundary1DVector(1, BoundaryElement(1.0));

  constexpr double inner_condition_value = -5;
  constexpr std::pair<bool, double> innerBoundary =
      std::make_pair(true, inner_condition_value);

  std::vector<std::pair<bool, double>> row_ibc(12, std::make_pair(false, -1));
  row_ibc[1] = innerBoundary;

  std::vector<std::pair<bool, double>> col_ibc(10, std::make_pair(false, -1));
  col_ibc[0] = innerBoundary;

  SUBCASE("Boundaries 1D case") {
    CHECK_NOTHROW(Boundary boundary(grid1D));
    CHECK_EQ(boundary1D.getBoundarySide(BC_SIDE_LEFT).size(), 1);
    CHECK_EQ(boundary1D.getBoundarySide(BC_SIDE_RIGHT).size(), 1);
    CHECK_EQ(boundary1D.getBoundaryElementType(BC_SIDE_LEFT, 0),
             BC_TYPE_CLOSED);
    CHECK_THROWS(boundary1D.getBoundarySide(BC_SIDE_TOP));
    CHECK_THROWS(boundary1D.getBoundarySide(BC_SIDE_BOTTOM));
    CHECK_NOTHROW(boundary1D.setBoundarySideClosed(BC_SIDE_LEFT));
    CHECK_THROWS(boundary1D.setBoundarySideClosed(BC_SIDE_TOP));
    CHECK_NOTHROW(boundary1D.setBoundarySideConstant(BC_SIDE_LEFT, 1.0));
    CHECK_EQ(boundary1D.getBoundaryElementValue(BC_SIDE_LEFT, 0), 1.0);
    CHECK_THROWS(boundary1D.getBoundaryElementValue(BC_SIDE_LEFT, 2));
    CHECK_EQ(boundary1D.getBoundaryElementType(BC_SIDE_LEFT, 0),
             BC_TYPE_CONSTANT);
    CHECK_EQ(boundary1D.getBoundaryElement(BC_SIDE_LEFT, 0).getType(),
             boundary1DVector[0].getType());

    CHECK_NOTHROW(boundary1D.setInnerBoundary(0, inner_condition_value));
    CHECK_THROWS(boundary1D.setInnerBoundary(0, 0, inner_condition_value));
    CHECK_EQ(boundary1D.getInnerBoundary(0), innerBoundary);
    CHECK_EQ(boundary1D.getInnerBoundary(1).first, false);
  }

  SUBCASE("Boundaries 2D case") {
    CHECK_NOTHROW(Boundary boundary(grid1D));
    CHECK_EQ(boundary2D.getBoundarySide(BC_SIDE_LEFT).size(), 10);
    CHECK_EQ(boundary2D.getBoundarySide(BC_SIDE_RIGHT).size(), 10);
    CHECK_EQ(boundary2D.getBoundarySide(BC_SIDE_TOP).size(), 12);
    CHECK_EQ(boundary2D.getBoundarySide(BC_SIDE_BOTTOM).size(), 12);
    CHECK_EQ(boundary2D.getBoundaryElementType(BC_SIDE_LEFT, 0),
             BC_TYPE_CLOSED);
    CHECK_NOTHROW(boundary2D.getBoundarySide(BC_SIDE_TOP));
    CHECK_NOTHROW(boundary2D.getBoundarySide(BC_SIDE_BOTTOM));
    CHECK_NOTHROW(boundary2D.setBoundarySideClosed(BC_SIDE_LEFT));
    CHECK_NOTHROW(boundary2D.setBoundarySideClosed(BC_SIDE_TOP));
    CHECK_NOTHROW(boundary2D.setBoundarySideConstant(BC_SIDE_LEFT, 1.0));
    CHECK_EQ(boundary2D.getBoundaryElementValue(BC_SIDE_LEFT, 0), 1.0);
    CHECK_THROWS(boundary2D.getBoundaryElementValue(BC_SIDE_LEFT, 12));
    CHECK_EQ(boundary2D.getBoundaryElementType(BC_SIDE_LEFT, 0),
             BC_TYPE_CONSTANT);
    CHECK_EQ(boundary2D.getBoundaryElement(BC_SIDE_LEFT, 0).getType(),
             boundary1DVector[0].getType());

    CHECK_THROWS(boundary2D.setInnerBoundary(0, inner_condition_value));
    CHECK_NOTHROW(boundary2D.setInnerBoundary(0, 1, inner_condition_value));
    CHECK_EQ(boundary2D.getInnerBoundary(0, 1), innerBoundary);
    CHECK_EQ(boundary2D.getInnerBoundary(0, 2).first, false);

    CHECK_EQ(boundary2D.getInnerBoundaryRow(0), row_ibc);
    CHECK_EQ(boundary2D.getInnerBoundaryCol(1), col_ibc);
  }
}
