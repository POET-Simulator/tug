#include <Eigen/Eigen>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <files.hpp>
#include <tug/Diffusion.hpp>

using namespace tug;

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

template <class T, tug::APPROACH app> int doWork(int ngrid) {

  constexpr T dt = 10;

  // create a grid
  std::string name;

  if constexpr (std::is_same_v<T, double>) {
    name = "DOUBLE";
  } else if constexpr (std::is_same_v<T, float>) {
    name = "FLOAT";
  } else {
    name = "unknown";
  }

  std::cout << name << " grid: " << ngrid << std::endl;

  Grid<T> grid(ngrid, ngrid);
  // Grid64 grid(ngrid, ngrid);

  // (optional) set the domain, e.g.:
  grid.setDomain(0.1, 0.1);

  Eigen::MatrixX<T> initConc64 = Eigen::MatrixX<T>::Constant(ngrid, ngrid, 0);
  initConc64(50, 50) = 1E-6;
  grid.setConcentrations(initConc64);

  const T sum_init64 = initConc64.sum();

  constexpr T alphax_val = 5e-10;
  Eigen::MatrixX<T> alphax =
      Eigen::MatrixX<T>::Constant(ngrid, ngrid, alphax_val); // row,col,value
  constexpr T alphay_val = 1e-10;
  Eigen::MatrixX<T> alphay =
      Eigen::MatrixX<T>::Constant(ngrid, ngrid, alphay_val); // row,col,value
  grid.setAlpha(alphax, alphay);

  // create a boundary with constant values
  Boundary bc = Boundary(grid);
  bc.setBoundarySideClosed(BC_SIDE_LEFT);
  bc.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc.setBoundarySideClosed(BC_SIDE_TOP);
  bc.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // set up a simulation environment
  Diffusion Sim(grid, bc); // grid_64,boundary,simulation-approach

  // Sim64.setSolver(THOMAS_ALGORITHM_SOLVER);

  // set the timestep of the simulation
  Sim.setTimestep(dt); // timestep

  // set the number of iterations
  Sim.setIterations(2);

  // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  Sim.setOutputCSV(CSV_OUTPUT_ON);

  // set output to the console to 'ON'
  Sim.setOutputConsole(CONSOLE_OUTPUT_OFF);

  // // **** RUN SIM64 ****
  auto begin_t = std::chrono::high_resolution_clock::now();

  Sim.run();

  auto end_t = std::chrono::high_resolution_clock::now();
  auto ms_t =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t);

  const double sum_after64 = grid.getConcentrations().sum();

  std::cout << "Sum of init field:      " << std::setprecision(15) << sum_init64
            << "\nSum after 2 iterations: " << sum_after64
            << "\nMilliseconds: " << ms_t.count() << std::endl
            << std::endl;
  return 0;
}

int main(int argc, char *argv[]) {

  int n[] = {101, 201, 501, 1001, 2001};

  for (int i = 0; i < std::size(n); i++) {
    doWork<float, tug::BTCS_APPROACH>(n[i]);
    doWork<double, tug::BTCS_APPROACH>(n[i]);
    doWork<float, tug::FTCS_APPROACH>(n[i]);
    doWork<double, tug::FTCS_APPROACH>(n[i]);
  }

  return 0;
}
