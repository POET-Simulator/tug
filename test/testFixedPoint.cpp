#include "tug/Datatypes/FixedPoint.hpp"

#include <cstdint>
#include <iostream>

using namespace tug;

using FPToTest = Fixed8<5>;

int main() {
  FPToTest with(1.5);
  FPToTest without(0.5);

  const int test = 0;

  const FPToTest lol = 0;

  if (lol == test) {
    std::cout << "This works\n";
  }

  double operation;

  operation = with + without;
  std::cout << "Addition: " << operation << std::endl;

  operation = with - without;
  std::cout << "Subtraction: " << operation << std::endl;

  operation = with * without;
  std::cout << "Mult: " << operation << std::endl;

  operation = with / without;
  std::cout << "Div: " << operation << std::endl;
}