#include <Eigen/Eigen>
#include <iostream>
#include <tug/Advection.hpp>
#include <tug/Core/Matrix.hpp>

using namespace Eigen;
using namespace tug;

int main(int argc, char *argv[]) {
  int row = 21;
  int col = 21;
  // create two grids of equal size, grid1 for hydraulics heads, grid2 for
  // Concentrations

  RowMajMat<double> hydHeads = RowMajMat<double>::Constant(row, col, 1);
  RowMajMat<double> concentrations = RowMajMat<double>::Constant(row, col, 0);

  Grid64 gridHeads(hydHeads);
  Grid64 gridConc(concentrations);
  gridHeads.setDomain(100, 100);
  gridConc.setDomain(100, 100);

  // create boundaries
  Boundary bcH = Boundary(gridHeads);
  bcH.setBoundarySideConstant(BC_SIDE_LEFT, 10);
  bcH.setBoundarySideConstant(BC_SIDE_RIGHT, 1);
  bcH.setBoundarySideClosed(BC_SIDE_TOP);
  bcH.setBoundarySideClosed(BC_SIDE_BOTTOM);
  Boundary bcC = Boundary(gridConc);
  bcC.setBoundarySideConstant(BC_SIDE_LEFT, 0.1);
  bcC.setBoundarySideConstant(BC_SIDE_RIGHT, 1);
  bcC.setBoundarySideClosed(BC_SIDE_TOP);
  bcC.setBoundarySideClosed(BC_SIDE_BOTTOM);

  Velocities velocities = Velocities(gridHeads, bcH);
  velocities.setOutputCSV(CSV_OUTPUT_ON);
  velocities.setK(1E-2);
  velocities.setEpsilon(1E-8);
  velocities.setInjh(10);
  velocities.setIterations(0);
  // calculate steady hydraulic heads
  velocities.run();

  std::cout << "Velocities simulation finished." << std::endl;
  std::cout << hydHeads << std::endl;

  // set true for steady case
  Advection advection = Advection(velocities, gridConc, bcC, true);
  advection.setPorosity(0.2);
  advection.setIterations(21);
  // set timestep close almost exactly to clf to test advection
  advection.setTimestep(5039.05);
  //   velocities.setOutputCSV(CSV_OUTPUT_VERBOSE);
  advection.run();

  std::cout << "Advection simulation finished." << std::endl;

  std::cout << concentrations << std::endl;
}