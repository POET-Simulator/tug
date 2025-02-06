#include "tug/Advection/Advection.hpp"
#include "tug/Advection/Velocities.hpp"
#include "tug/Core/Matrix.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <tug/tug.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {
  int row = 5;
  int col = 5;
  // create two grids of equal size, grid1 for hydraulics heads, grid2 for
  // Concentrations

  RowMajMat<double> hydHeads = RowMajMat<double>::Constant(row, col, 1);
  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);
  hydHeads(row / 2, col / 2) = 10;
  concentrations(row / 2, col / 2) = 1;

  Velocities<double, tug::HYDRAULIC_MODE::STEADY_STATE,
             tug::HYDRAULIC_RESOLVE::EXPLICIT>
      velocities(hydHeads);
  velocities.setDomain(100, 100);
  velocities.setPermKX(RowMajMat<double>::Constant(row, col, 1E-8));
  velocities.setPermKY(RowMajMat<double>::Constant(row, col, 1E-8));
  velocities.setEpsilon(1E-8);

  Advection advection = Advection(concentrations, velocities);

  advection.setPorosity(RowMajMat<double>::Constant(row, col, 0.2));
  advection.setIterations(1);
  advection.setTimestep(10000);

  // create boundaries
  Boundary<double> &bcH = velocities.getBoundaryConditions();
  bcH.setBoundarySideConstant(BC_SIDE_LEFT, 1);
  bcH.setBoundarySideConstant(BC_SIDE_RIGHT, 1);
  bcH.setBoundarySideConstant(BC_SIDE_TOP, 1);
  bcH.setBoundarySideConstant(BC_SIDE_BOTTOM, 1);
  // bcH.setBoundarySideClosed(BC_SIDE_TOP);
  // bcH.setBoundarySideClosed(BC_SIDE_BOTTOM);
  bcH.setInnerBoundary(row / 2, col / 2, 10);
  //
  Boundary<double> &bcC = advection.getBoundaryConditions();
  // bcC.setBoundarySideConstant(BC_SIDE_LEFT, 0.1);
  // bcC.setBoundarySideConstant(BC_SIDE_RIGHT, 1);
  bcC.setBoundarySideClosed(BC_SIDE_TOP);
  bcC.setBoundarySideClosed(BC_SIDE_BOTTOM);

  advection.run();

  std::cout << velocities.getConcentrationMatrix() << std::endl;
  std::cout << velocities.getVelocitiesX() << std::endl;
  std::cout << velocities.getVelocitiesY() << std::endl;

  std::cout << "Advection simulation finished." << std::endl;

  std::cout << concentrations << std::endl;
}
