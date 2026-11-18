#include <cstddef>
#include <iostream>
#include <tug/tug.hpp>

using namespace tug;

int main(int argc, char *argv[]) {
  constexpr std::size_t row = 5;
  constexpr std::size_t col = 5;

  RowMajMat<double> hydHeads = RowMajMat<double>::Constant(row, col, 1);
  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);

  Velocities<double, HYDRAULIC_MODE::STEADY_STATE, HYDRAULIC_RESOLVE::EXPLICIT>
      velocities(hydHeads);

  velocities.setDomain(5, 5);
  velocities.setPermKX(RowMajMat<double>::Constant(row, col, 3E-7));
  velocities.setPermKY(RowMajMat<double>::Constant(row, col, 3E-7));
  velocities.setEpsilon(1E-8);

  Advection advection(concentrations, velocities);

  advection.setPorosity(0.2);
  advection.setIterations(3);
  advection.setTimestep(3.3333e+05);

  // Set hydraulic boundaries
  Boundary<double> &bcH = velocities.getBoundaryConditions();
  bcH.setBoundarySideConstant(BC_SIDE_LEFT, 10);
  bcH.setBoundarySideConstant(BC_SIDE_RIGHT, 0);

  // Set concentration boundaries
  Boundary<double> &bcC = advection.getBoundaryConditions();
  bcC.setBoundarySideConstant(BC_SIDE_LEFT, 1);
  bcC.setBoundarySideConstant(BC_SIDE_RIGHT, 0);

  advection.run();

  std::cout << "Hydraulic Heads:\n"
            << velocities.getConcentrationMatrix() << "\n\n";

  std::cout << "Velocities X:\n"
            << velocities.getVelocitiesX() << "\n\nVelocities Y:\n"
            << velocities.getVelocitiesY() << "\n\n";

  std::cout << "Concentrations:\n" << concentrations << std::endl;

  return 0;
}
