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

/**
 * Try to parse an input string into a given template type.
 */
template <typename T> inline T parseString(const std::string &str) {
  T result;
  std::istringstream iss(str);

  if (!(iss >> result)) {
    throw std::invalid_argument("Invalid input for parsing.");
  }

  return result;
}

/**
 * Splits a given string into a vector by using a delimiter character.
 */
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

/**
 * Opens a file containing CSV and transform it into row-major 2D STL vector.
 */
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

/**
 * Converts a 2D STL vector, where values are stored row-major into a
 * column-major Eigen::Matrix.
 */
template <typename T>
Eigen::MatrixXd rmVecTocmMatrix(const std::vector<std::vector<T>> &vec,
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
  Grid64 grid(row, col);

  // (optional) set the domain, e.g.:
  grid.setDomain(0.005, 0.01);

  const auto init_values_vec = CSVToVector<double>(INPUT_CONC_FILE);
  Eigen::MatrixXd concentrations = rmVecTocmMatrix(init_values_vec, row, col);
  grid.setConcentrations(concentrations);

  const double sum_init = concentrations.sum();

  // // (optional) set alphas of the grid, e.g.:
  const auto alphax_vec = CSVToVector<double>(INPUT_ALPHAX_FILE);
  Eigen::MatrixXd alphax = rmVecTocmMatrix(alphax_vec, row, col);

  constexpr double alphay_val = 5e-10;
  Eigen::MatrixXd alphay = Eigen::MatrixXd::Constant(row, col, alphay_val); // row,col,value
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
  Simulation simulation =
      Simulation(grid, bc, BTCS_APPROACH); // grid,boundary,simulation-approach

  // set the timestep of the simulation
  simulation.setTimestep(360); // timestep

  // set the number of iterations
  simulation.setIterations(1);

  // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  simulation.setOutputCSV(CSV_OUTPUT_ON);

  // set output to the console to 'ON'
  simulation.setOutputConsole(CONSOLE_OUTPUT_ON);

  // // **** RUN SIMULATION ****
  simulation.run();

  const double sum_after = grid.getConcentrations().sum();

  std::cout << "Sum of init field: " << std::to_string(sum_init)
            << "\nSum after iteration: " << std::to_string(sum_after)
            << std::endl;

  return 0;
}
