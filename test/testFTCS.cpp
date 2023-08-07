#include <doctest/doctest.h>
#include <../src/FTCS.cpp>
#include <limits>

TEST_CASE("Maths") {
    SUBCASE("mean between two alphas") {
        double alpha1 = 10;
        double alpha2 = 20;
        double average = 15;
        double harmonicMean = double(2) / ((double(1)/alpha1)+(double(1)/alpha2));

        // double difference = std::fabs(calcAlphaIntercell(alpha1, alpha2) - harmonicMean);
        // CHECK(difference < std::numeric_limits<double>::epsilon());
        CHECK_EQ(calcAlphaIntercell(alpha1, alpha2), harmonicMean);
        CHECK_EQ(calcAlphaIntercell(alpha1, alpha2, false), average);
    }


}