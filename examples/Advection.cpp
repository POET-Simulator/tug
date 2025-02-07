#include "tug/Boundary.hpp"
#include <cstddef>
#include <iostream>
#include <tug/tug.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {
  constexpr std::size_t row = 5;
  constexpr std::size_t col = 5;

  RowMajMat<double> hydHeads = RowMajMat<double>::Constant(row, col, 1);
  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);

  Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
             tug::HYDRAULIC_RESOLVE::EXPLICIT>
      velocities(hydHeads);

  velocities.setDomain(1, 1);
  velocities.setPermKX(RowMajMat<double>::Constant(row, col, 3E-7));
  velocities.setPermKY(RowMajMat<double>::Constant(row, col, 3E-7));
  velocities.setEpsilon(1E-8);

  Advection advection = Advection(concentrations, velocities);

  advection.setPorosity(RowMajMat<double>::Constant(row, col, 0.2));
  advection.setIterations(6);
  // 1 hour
  advection.setTimestep(6666);

  // create boundaries
  Boundary<double> &bcH = velocities.getBoundaryConditions();
  bcH.setBoundarySideConstant(BC_SIDE_LEFT, 10);
  bcH.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
  // bcH.setBoundarySideConstant(BC_SIDE_TOP, 1);
  // bcH.setBoundarySideConstant(BC_SIDE_BOTTOM, 1);
  // bcH.setInnerBoundary(row / 2, col / 2, 10);

  Boundary<double> &bcC = advection.getBoundaryConditions();
  bcC.setBoundarySideConstant(BC_SIDE_LEFT, 1);
  bcC.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
  // bcC.setInnerBoundary(row / 2, col / 2, 1);
  // bcC.setBoundarySideConstant(BC_SIDE_LEFT, 0);
  // bcC.setBoundarySideConstant(BC_SIDE_RIGHT, 0);
  // bcC.setBoundarySideConstant(BC_SIDE_TOP, 0);
  // bcC.setBoundarySideConstant(BC_SIDE_BOTTOM, 0);

  advection.run();

  std::cout << velocities.getConcentrationMatrix() << std::endl << std::endl;

  std::cout << velocities.getVelocitiesX() << std::endl
            << std::endl
            << velocities.getVelocitiesY() << std::endl
            << std::endl;

  std::cout << concentrations << std::endl;
}
