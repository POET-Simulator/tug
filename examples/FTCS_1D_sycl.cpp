#include <chrono>
#include <tug/Grid.hpp>
#include <tug/Simulation.hpp>

#include <Eigen/Dense>

using namespace tug;
using namespace Eigen;

MatrixXd runNormal(MatrixXd concentrations, const MatrixXd &alphas,
                   Boundary<double> &boundaries, int iterations,
                   double timestep) {
  Grid64 grid(concentrations.cols());
  grid.setConcentrations(concentrations);

  grid.setAlpha(alphas);

  Simulation simulation =
      Simulation<double, tug::FTCS_APPROACH>(grid, boundaries);

  simulation.setTimestep(timestep);
  simulation.setIterations(iterations);
  simulation.run();

  return MatrixXd(grid.getConcentrations());
}

MatrixXd runWithSYCL(MatrixXd concentrations, const MatrixXd &alphas,
                     Boundary<double> &boundaries, int iterations,
                     double timestep) {
  Grid64 grid(concentrations.cols());
  grid.setConcentrations(concentrations);

  grid.setAlpha(alphas);

  Simulation simulation = Simulation<double, tug::FTCS_SYCL>(grid, boundaries);

  simulation.setTimestep(timestep);
  simulation.setIterations(iterations);
  simulation.run();

  return MatrixXd(grid.getConcentrations());
}

int main(int argc, char **argv) {
  constexpr int cells = 1e7;
  constexpr double timestep = 0.1;
  constexpr int iterations = 500;

  Grid64 grid(cells);

  MatrixXd concentrations = MatrixXd::Constant(1, cells, 20);

  MatrixXd alphas = MatrixXd::Constant(1, cells, 1);

  // ******************
  // **** BOUNDARY ****
  // ******************

  // create a boundary with constant values
  Boundary bc = Boundary(grid);
  bc.setBoundarySideConstant(BC_SIDE_LEFT, 1);
  bc.setBoundarySideConstant(BC_SIDE_RIGHT, 1);

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

  auto diff = (normal - sycl).cwiseAbs().maxCoeff();

  std::cout << "Max diff: " << diff << std::endl;

  return 0;
}