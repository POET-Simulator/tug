#include <Eigen/Eigen>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tug/Simulation.hpp>
#include <vector>

#include "files.hpp"

template <typename T> inline T parseString(const std::string &str) {
  T result;
  std::istringstream iss(str);

  if (!(iss >> result)) {
    throw std::invalid_argument("Invalid input for parsing.");
  }

  return result;
}

template <typename T>
std::vector<T> tokenize(const std::string &input, char delimiter) {
  std::vector<T> tokens;
  std::istringstream tokenStream(input);
  std::string token;

  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(parseString<T>(token));
  }

  return tokens;
}

template <typename T>
std::vector<std::vector<T>> CSVToVector(const char *filename) {
  std::ifstream in_file(filename);
  if (!in_file.is_open()) {
    throw std::runtime_error("Error opening file \'" + std::string(filename) +
                             "\'.");
  }

  std::vector<std::vector<T>> csv_data;

  std::string line;
  while (std::getline(in_file, line)) {
    csv_data.push_back(tokenize<T>(line, ','));
  }

  in_file.close();

  return csv_data;
}

template <typename T>
Eigen::MatrixXd CMVecToRMMatrix(const std::vector<std::vector<T>> &vec,
                                std::uint32_t exp_rows,
                                std::uint32_t exp_cols) {
  if (exp_rows != vec.size()) {
    throw std::runtime_error(
        "Mismatch in y dimension while converting to Eigen::Matrix.");
  }

  Eigen::MatrixXd out_mat(exp_rows, exp_cols);

  for (std::uint32_t ri = 0; ri < exp_rows; ri++) {
    const auto &vec_row = vec[ri];
    if (vec[ri].size() != exp_cols) {
      throw std::runtime_error(
          "Mismatch in x dimension while converting to Eigen::Matrix.");
    }
    for (std::uint32_t cj = 0; cj < exp_cols; cj++) {
      out_mat(ri, cj) = vec_row[cj];
    }
  }

  return out_mat;
}

int main(int argc, char *argv[]) {
  // EASY_PROFILER_ENABLE;
  // profiler::startListen();
  // **************
  // **** GRID ****
  // **************
  // profiler::startListen();
  // create a grid with a 5 x 10 field
  constexpr int row = 5;
  constexpr int col = 10;
  Grid grid = Grid(row, col);

  // (optional) set the domain, e.g.:
  grid.setDomain(0.005, 0.01);

  const auto init_values_vec = CSVToVector<double>(INPUT_CONC_FILE);
  MatrixXd concentrations = CMVecToRMMatrix(init_values_vec, row, col);
  grid.setConcentrations(concentrations);

  // // (optional) set alphas of the grid, e.g.:
  const auto alphax_vec = CSVToVector<double>(INPUT_ALPHAX_FILE);
  MatrixXd alphax = CMVecToRMMatrix(alphax_vec, row, col);

  constexpr double alphay_val = 5e-10;
  MatrixXd alphay = MatrixXd::Constant(row, col, alphay_val); // row,col,value
  grid.setAlpha(alphax, alphay);

  // // ******************
  // // **** BOUNDARY ****
  // // ******************

  // create a boundary with constant values
  Boundary bc = Boundary(grid);
  bc.setBoundarySideClosed(BC_SIDE_LEFT);
  bc.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc.setBoundarySideClosed(BC_SIDE_TOP);
  bc.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // // ************************
  // // **** SIMULATION ENV ****
  // // ************************

  // set up a simulation environment
  Simulation simulation = Simulation(
      grid, bc, BTCS_APPROACH); // grid,boundary,simulation-approach

  // set the timestep of the simulation
  simulation.setTimestep(360); // timestep

  // set the number of iterations
  simulation.setIterations(1);

  // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT_ON);

  simulation.setOutputConsole(CONSOLE_OUTPUT_ON);

  // // **** RUN SIMULATION ****

  // // run the simulation

  // // EASY_BLOCK("SIMULATION")
  simulation.run();
  // // EASY_END_BLOCK;
  // // profiler::dumpBlocksToFile("test_profile.prof");
  // // profiler::stopListen();
}
