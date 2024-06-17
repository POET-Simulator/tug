#include <chrono>
#include <tug/Grid.hpp>
#include <tug/Simulation.hpp>

#include <Eigen/Dense>

using namespace tug;
using namespace Eigen;

MatrixXf runNormal(MatrixXf concentrations, const MatrixXf &alphas,
                   Boundary<float> &boundaries, int iterations,
                   double timestep) {
  Grid32 grid(concentrations.rows(), concentrations.cols());
  grid.setConcentrations(concentrations);

  grid.setAlpha(alphas, alphas);

  Simulation simulation =
      Simulation<float, tug::FTCS_APPROACH>(grid, boundaries);

  simulation.setTimestep(timestep);
  simulation.setIterations(iterations);
  simulation.run();

  return MatrixXf(grid.getConcentrations());
}

MatrixXf runWithSYCL(MatrixXf concentrations, const MatrixXf &alphas,
                     Boundary<float> &boundaries, int iterations,
                     double timestep) {
  Grid32 grid(concentrations.rows(), concentrations.cols());
  grid.setConcentrations(concentrations);

  grid.setAlpha(alphas, alphas);

  Simulation simulation = Simulation<float, tug::FTCS_SYCL>(grid, boundaries);

  simulation.setTimestep(timestep);
  simulation.setIterations(iterations);
  simulation.run();

  return MatrixXf(grid.getConcentrations());
}

int main(int argc, char **argv) {
  constexpr int cells = 3e3;
  constexpr double timestep = 500;
  constexpr int iterations = 1;

  MatrixXf concentrations(cells, cells);
  concentrations.setConstant(20);
  concentrations(0, 0) = 100;

  MatrixXf alphas(cells, cells);
  alphas.setConstant(1);

  // ******************
  // **** BOUNDARY ****
  // ******************

  // create a boundary with constant values
  Boundary bc = Boundary<float>(cells, cells);
  // bc.setBoundarySideConstant(BC_SIDE_LEFT, 1);
  // bc.setBoundarySideConstant(BC_SIDE_RIGHT, 1);

  auto t_normal_begin = std::chrono::steady_clock::now();
  auto normal = runNormal(concentrations, alphas, bc, iterations, timestep);
  auto t_normal_end = std::chrono::steady_clock::now();

  std::cout << "Normal time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   t_normal_end - t_normal_begin)
                   .count()
            << "ms" << std::endl;

  auto t_sycl_begin = std::chrono::steady_clock::now();
  auto sycl = runWithSYCL(concentrations, alphas, bc, iterations, timestep);
  auto t_sycl_end = std::chrono::steady_clock::now();

  std::cout << "SYCL time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(
                   t_sycl_end - t_sycl_begin)
                   .count()
            << "ms" << std::endl;

  // std::cout << "Normal: " << std::endl << normal << std::endl;
  // std::cout << "SYCL: " << std::endl << sycl << std::endl;

  auto diff = normal - sycl;

  std::cout << "Max diff: " << diff.cwiseAbs().maxCoeff() << std::endl;

  return 0;
}