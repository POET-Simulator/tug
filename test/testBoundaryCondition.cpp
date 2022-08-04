#include <grid/BTCSBoundaryCondition.hpp>
#include <doctest/doctest.h>

using namespace Diffusion;

#define BC_CONST_VALUE 1e-5

TEST_CASE("1D Boundary Condition") {

  BTCSBoundaryCondition bc(5);
  boundary_condition bc_set = {BC_TYPE_CONSTANT, BC_CONST_VALUE};

  SUBCASE("valid get") { CHECK_EQ(bc(BC_SIDE_LEFT).value, 0); }

  SUBCASE("invalid get") {
    CHECK_THROWS(bc(BC_SIDE_TOP));
    CHECK_THROWS(bc(BC_SIDE_LEFT, 1));
  }

  SUBCASE("valid set") {
    CHECK_NOTHROW(bc(BC_SIDE_LEFT) = bc_set);
    CHECK_EQ(bc(BC_SIDE_LEFT).value, bc_set.value);
    CHECK_EQ(bc(BC_SIDE_LEFT).type, bc_set.type);
  }

  SUBCASE("invalid set") {
    CHECK_THROWS(bc(BC_SIDE_TOP) = bc_set);
    CHECK_THROWS(bc(BC_SIDE_LEFT, 0) = bc_set);
  }

  SUBCASE("valid row getter") {
    bc(BC_SIDE_LEFT) = bc_set;
    bc_tuple tup = bc.row_boundary(0);

    CHECK_EQ(tup[0].value, BC_CONST_VALUE);
    CHECK_EQ(tup[1].value, 0);
  }

  SUBCASE("invalid row and col getter") {
    CHECK_THROWS(bc.row_boundary(1));
    CHECK_THROWS(bc.col_boundary(0));
  }
}

TEST_CASE("2D Boundary Condition") {

  BTCSBoundaryCondition bc(5, 5);
  boundary_condition bc_set = {BC_TYPE_CONSTANT, BC_CONST_VALUE};

  SUBCASE("valid get") { CHECK_EQ(bc(BC_SIDE_LEFT, 0).value, 0); }

  SUBCASE("invalid get") {
    CHECK_THROWS(bc(5, 0));
    CHECK_THROWS(bc(BC_SIDE_LEFT));
  }

  SUBCASE("valid set") {
    CHECK_NOTHROW(bc(BC_SIDE_LEFT, 0) = bc_set);
    CHECK_EQ(bc(BC_SIDE_LEFT, 0).value, bc_set.value);
    CHECK_EQ(bc(BC_SIDE_LEFT, 0).type, bc_set.type);
  }

  SUBCASE("invalid set") {
    CHECK_THROWS(bc(BC_SIDE_LEFT) = bc_set);
    CHECK_THROWS(bc(5, 0) = bc_set);
  }

  SUBCASE("call of setSide") {
    CHECK_NOTHROW(bc.setSide(BC_SIDE_BOTTOM, bc_set));
    CHECK_EQ(bc(BC_SIDE_BOTTOM, 1).value, bc_set.value);
    CHECK_EQ(bc(BC_SIDE_BOTTOM, 1).type, bc_set.type);
  }

  SUBCASE("get and set of side") {
    std::vector<boundary_condition> bc_vec;
    CHECK_NOTHROW(bc_vec = bc.getSide(BC_SIDE_BOTTOM));
    bc_vec[3] = {BC_TYPE_CONSTANT, 1e-5};
    CHECK_NOTHROW(bc.setSide(BC_SIDE_BOTTOM, bc_vec));
    CHECK_EQ(bc(BC_SIDE_BOTTOM, 3).type, BC_TYPE_CONSTANT);
    CHECK_EQ(bc(BC_SIDE_BOTTOM, 3).value, 1e-5);

    CHECK_EQ(bc(BC_SIDE_BOTTOM, 2).value, 0);
  }
}

TEST_CASE("Boundary Condition helpers") {
  boundary_condition bc_set = {BC_TYPE_CONSTANT, BC_CONST_VALUE};

  SUBCASE("return boundary condition skeleton") {
    boundary_condition bc_test = BTCSBoundaryCondition::returnBoundaryCondition(
        bc_set.type, bc_set.value);
    CHECK_EQ(bc_test.value, bc_set.value);
    CHECK_EQ(bc_test.type, bc_set.type);
  }
}
