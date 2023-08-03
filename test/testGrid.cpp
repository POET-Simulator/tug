#include <doctest/doctest.h>
#include <tug/Grid.hpp>

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

    SUBCASE("") {

    }
}

TEST_CASE("2D Grid quadratic") {
    int r = 12;
    int c = 12;

    // TODO
}

TEST_CASE("2D Grid non-quadratic") {
    int r = 12;
    int c = 15;

    // TODO
}