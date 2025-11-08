#include "bench_defs.hpp"
#include "io.hpp"
#include <Boundary.hpp>
#include <Diffusion.hpp>
#include <bench_defs.hpp>
#include <cstdlib>
#include <iostream>

using TugType = double;

void run_bench(const bench_input &input, const std::string &output_file) {

  std::vector<std::vector<TugType>> raw_data =
      read_conc_csv<TugType>(input.csv_file_init, input.nrows, input.ncols);
  std::vector<TugType> alpha_x_input =
      read_alpha_csv<TugType>(input.csv_alpha_x);
  std::vector<TugType> alpha_y_input =
      read_alpha_csv<TugType>(input.csv_alpha_y);

  Eigen::Map<tug::RowMajMat<TugType>> alpha_x(alpha_x_input.data(), input.nrows,
                                              input.ncols);
  Eigen::Map<tug::RowMajMat<TugType>> alpha_y(alpha_y_input.data(), input.nrows,
                                              input.ncols);

  for (int i = 0; i < raw_data.size(); i++) {

    tug::Diffusion<TugType> diffusion(raw_data[i].data(), input.nrows,
                                      input.ncols);
    diffusion.setAlphaX(alpha_x);
    diffusion.setAlphaY(alpha_y);
    diffusion.setDomain(input.s_x, input.s_y);
    diffusion.setIterations(input.iterations);
    diffusion.setTimestep(input.timestep);

    tug::Boundary<TugType> &bc = diffusion.getBoundaryConditions();

    // set north boundary
    for (const auto &index : input.boundary.north_const) {
      bc.setBoundaryElementConstant(tug::BC_SIDE_TOP, index,
                                    input.boundary.values[i]);
    }

    // set south boundary
    for (const auto &index : input.boundary.south_const) {
      bc.setBoundaryElementConstant(tug::BC_SIDE_BOTTOM, index,
                                    input.boundary.values[i]);
    }

    // set east boundary
    for (const auto &index : input.boundary.east_const) {
      bc.setBoundaryElementConstant(tug::BC_SIDE_RIGHT, index,
                                    input.boundary.values[i]);
    }

    // set west boundary
    for (const auto &index : input.boundary.west_const) {
      bc.setBoundaryElementConstant(tug::BC_SIDE_LEFT, index,
                                    input.boundary.values[i]);
    }

    diffusion.run();
  }

  write_conc_csv(output_file, raw_data);
  std::cout << "Output written to: " << output_file << std::endl;
}

int main(int argc, char *argv[]) {
  int input_index = 0;
  if (argc > 1) {
    input_index = std::atoi(argv[1]);
  } else {
    std::cout << "Usage: " << argv[0] << " [input_index]" << std::endl;
    std::cout << "Available inputs:" << std::endl;
    std::cout << "  0 - barite_200_input (default)" << std::endl;
    std::cout << "  1 - barite_large_input" << std::endl;
    std::cout << "  2 - surfex_input" << std::endl;
  }

  switch (input_index) {
  case 0:
    std::cout << "run barite 200" << std::endl;
    run_bench(barite_200_input, "result_bartite_200.csv");
    break;
  case 1:
    std::cout << "run barite large" << std::endl;
    run_bench(barite_large_input, "result_barite_large.csv");
    break;
  case 2:
    std::cout << "run surfex" << std::endl;
    run_bench(surfex_input, "result_surfex.csv");
    break;
  default:
    std::cout << "Invalid input index " << input_index
              << ", using default (barite_200_input)" << std::endl;
    run_bench(barite_200_input, "result_barite_200.csv");
    break;
  }
  return 0;
}