#include <tug/Boundary.hpp>
#include <utility>
#include <vector>

using namespace std;
using namespace tug;

#include <gtest/gtest.h>

#define BOUNDARY_TEST(x) TEST(Boundary, x)

BOUNDARY_TEST(Element) {

  BoundaryElement boundaryElementClosed = BoundaryElement<double>();
  EXPECT_NO_THROW(BoundaryElement<double>());
  EXPECT_EQ(boundaryElementClosed.getType(), BC_TYPE_CLOSED);
  EXPECT_DOUBLE_EQ(boundaryElementClosed.getValue(), -1);
  EXPECT_ANY_THROW(boundaryElementClosed.setValue(0.2));

  BoundaryElement boundaryElementConstant = BoundaryElement(0.1);
  EXPECT_NO_THROW(BoundaryElement(0.1));
  EXPECT_EQ(boundaryElementConstant.getType(), BC_TYPE_CONSTANT);
  EXPECT_DOUBLE_EQ(boundaryElementConstant.getValue(), 0.1);
  EXPECT_NO_THROW(boundaryElementConstant.setValue(0.2));
  EXPECT_DOUBLE_EQ(boundaryElementConstant.getValue(), 0.2);
}

BOUNDARY_TEST(Class) {
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

  {
    EXPECT_NO_THROW(Boundary boundary(grid1D));
    EXPECT_EQ(boundary1D.getBoundarySide(BC_SIDE_LEFT).size(), 1);
    EXPECT_EQ(boundary1D.getBoundarySide(BC_SIDE_RIGHT).size(), 1);
    EXPECT_EQ(boundary1D.getBoundaryElementType(BC_SIDE_LEFT, 0),
              BC_TYPE_CLOSED);
    EXPECT_ANY_THROW(boundary1D.getBoundarySide(BC_SIDE_TOP));
    EXPECT_ANY_THROW(boundary1D.getBoundarySide(BC_SIDE_BOTTOM));
    EXPECT_NO_THROW(boundary1D.setBoundarySideClosed(BC_SIDE_LEFT));
    EXPECT_ANY_THROW(boundary1D.setBoundarySideClosed(BC_SIDE_TOP));
    EXPECT_NO_THROW(boundary1D.setBoundarySideConstant(BC_SIDE_LEFT, 1.0));
    EXPECT_DOUBLE_EQ(boundary1D.getBoundaryElementValue(BC_SIDE_LEFT, 0), 1.0);
    EXPECT_ANY_THROW(boundary1D.getBoundaryElementValue(BC_SIDE_LEFT, 2));
    EXPECT_EQ(boundary1D.getBoundaryElementType(BC_SIDE_LEFT, 0),
              BC_TYPE_CONSTANT);
    EXPECT_EQ(boundary1D.getBoundaryElement(BC_SIDE_LEFT, 0).getType(),
              boundary1DVector[0].getType());

    EXPECT_NO_THROW(boundary1D.setInnerBoundary(0, inner_condition_value));
    EXPECT_ANY_THROW(boundary1D.setInnerBoundary(0, 0, inner_condition_value));
    EXPECT_EQ(boundary1D.getInnerBoundary(0), innerBoundary);
    EXPECT_FALSE(boundary1D.getInnerBoundary(1).first);
  }

  {
    EXPECT_NO_THROW(Boundary boundary(grid1D));
    EXPECT_EQ(boundary2D.getBoundarySide(BC_SIDE_LEFT).size(), 10);
    EXPECT_EQ(boundary2D.getBoundarySide(BC_SIDE_RIGHT).size(), 10);
    EXPECT_EQ(boundary2D.getBoundarySide(BC_SIDE_TOP).size(), 12);
    EXPECT_EQ(boundary2D.getBoundarySide(BC_SIDE_BOTTOM).size(), 12);
    EXPECT_EQ(boundary2D.getBoundaryElementType(BC_SIDE_LEFT, 0),
              BC_TYPE_CLOSED);
    EXPECT_NO_THROW(boundary2D.getBoundarySide(BC_SIDE_TOP));
    EXPECT_NO_THROW(boundary2D.getBoundarySide(BC_SIDE_BOTTOM));
    EXPECT_NO_THROW(boundary2D.setBoundarySideClosed(BC_SIDE_LEFT));
    EXPECT_NO_THROW(boundary2D.setBoundarySideClosed(BC_SIDE_TOP));
    EXPECT_NO_THROW(boundary2D.setBoundarySideConstant(BC_SIDE_LEFT, 1.0));
    EXPECT_DOUBLE_EQ(boundary2D.getBoundaryElementValue(BC_SIDE_LEFT, 0), 1.0);
    EXPECT_ANY_THROW(boundary2D.getBoundaryElementValue(BC_SIDE_LEFT, 12));
    EXPECT_EQ(boundary2D.getBoundaryElementType(BC_SIDE_LEFT, 0),
              BC_TYPE_CONSTANT);
    EXPECT_EQ(boundary2D.getBoundaryElement(BC_SIDE_LEFT, 0).getType(),
              boundary1DVector[0].getType());

    EXPECT_ANY_THROW(boundary2D.setInnerBoundary(0, inner_condition_value));
    EXPECT_NO_THROW(boundary2D.setInnerBoundary(0, 1, inner_condition_value));
    EXPECT_EQ(boundary2D.getInnerBoundary(0, 1), innerBoundary);
    EXPECT_FALSE(boundary2D.getInnerBoundary(0, 2).first);

    EXPECT_EQ(boundary2D.getInnerBoundaryRow(0), row_ibc);
    EXPECT_EQ(boundary2D.getInnerBoundaryCol(1), col_ibc);
  }
}
