#include <grid/BTCSBoundaryCondition.hpp>
#include <doctest/doctest.h>

using namespace tug::boundary_condition;

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

TEST_CASE("1D special inner grid cells") {
  BTCSBoundaryCondition bc(5);
  boundary_condition bc_set = {BC_TYPE_CONSTANT, BC_CONST_VALUE};

  SUBCASE("valid set") {
    CHECK_NOTHROW(bc(BC_INNER, 0) = bc_set);
  }

  SUBCASE("valid get") {
    bc(BC_INNER, 0) = bc_set;
    CHECK_EQ(bc(BC_INNER, 0).type, bc_set.type);
  }

  SUBCASE("invalid get") {
    CHECK_EQ(bc(BC_INNER, 1).type, BC_UNSET);
    CHECK_THROWS(bc(BC_INNER, 5));
  }

  SUBCASE("invalid set") {
    CHECK_THROWS(bc(BC_INNER, 5) = bc_set);
  }

  SUBCASE("valid row getter") {
    bc(BC_INNER, 1) = bc_set;

    bc_vec ret = bc.getInnerRow(0);

    CHECK_EQ(ret[0].type, BC_UNSET);
    CHECK_EQ(ret[1].type, bc_set.type);
  }

  SUBCASE("invalid row getter") {
    CHECK_THROWS(bc.getInnerRow(1));
  }

  SUBCASE("invalid col getter") {
    CHECK_THROWS(bc.getInnerCol(0));
  }

}

TEST_CASE("2D special inner grid cells") {
  BTCSBoundaryCondition bc(5,5);
  boundary_condition bc_set = {BC_TYPE_CONSTANT, BC_CONST_VALUE};

  SUBCASE("valid set") {
    CHECK_NOTHROW(bc(BC_INNER, 0) = bc_set);
  }

  SUBCASE("valid get") {
    bc(BC_INNER, 0) = bc_set;
    CHECK_EQ(bc(BC_INNER, 0).type, bc_set.type);
  }

  SUBCASE("invalid get") {
    CHECK_EQ(bc(BC_INNER, 1).type, BC_UNSET);
    CHECK_THROWS(bc(BC_INNER, 25));
  }

  SUBCASE("invalid set") {
    CHECK_THROWS(bc(BC_INNER, 25) = bc_set);
  }

  SUBCASE("valid row getter") {
    bc(BC_INNER, 0) = bc_set;
    bc(BC_INNER, 6) = bc_set;

    bc_vec ret = bc.getInnerRow(0);

    CHECK_EQ(ret[0].type, bc_set.type);
    CHECK_EQ(ret[1].type, BC_UNSET);

    ret = bc.getInnerRow(1);

    CHECK_EQ(ret[0].type, BC_UNSET);
    CHECK_EQ(ret[1].type, bc_set.type);
  }

  SUBCASE("valid col getter") {
    bc(BC_INNER, 1) = bc_set;
    bc(BC_INNER, 5) = bc_set;

    bc_vec ret = bc.getInnerCol(0);

    CHECK_EQ(ret[0].type, BC_UNSET);
    CHECK_EQ(ret[1].type, bc_set.type);

    ret = bc.getInnerCol(1);

    CHECK_EQ(ret[0].type, bc_set.type);
    CHECK_EQ(ret[1].type, BC_UNSET);
  }

  SUBCASE("invalid row getter") {
    CHECK_THROWS(bc.getInnerRow(5));
  }

  SUBCASE("invalid col getter") {
    CHECK_THROWS(bc.getInnerCol(5));
  }
}
