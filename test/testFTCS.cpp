#include <gtest/gtest.h>
#include <tug/Core/TugUtils.hpp>

#include <gtest/gtest.h>

TEST(FTCS, calcAlphaIntercell) {
  double alpha1 = 10;
  double alpha2 = 20;
  double average = 15;
  double harmonicMean =
      double(2) / ((double(1) / alpha1) + (double(1) / alpha2));

  // double difference = std::fabs(calcAlphaIntercell(alpha1, alpha2) -
  // harmonicMean); CHECK(difference <
  // std::numeric_limits<double>::epsilon());
  EXPECT_DOUBLE_EQ(calcAlphaIntercell(alpha1, alpha2), harmonicMean);
  EXPECT_DOUBLE_EQ(calcAlphaIntercell(alpha1, alpha2, false), average);
}
