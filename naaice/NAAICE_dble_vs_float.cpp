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
#include <chrono>

#include "files.hpp"

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


int DoDble(int ngrid, APPROACH approach) {
  
  constexpr double dt = 10;
  
  // create a grid
  std::cout << "DOUBLE grid: " << ngrid << std::endl;

  Grid64 grid64(ngrid, ngrid);

  // (optional) set the domain, e.g.:
  grid64.setDomain(0.1, 0.1);

  Eigen::MatrixXd initConc64 = Eigen::MatrixXd::Constant(ngrid, ngrid, 0);
  initConc64(50, 50) = 1E-6;
  grid64.setConcentrations(initConc64);
  
  const double sum_init64 = initConc64.sum();

  constexpr double alphax_val = 5e-10;
  Eigen::MatrixXd alphax = Eigen::MatrixXd::Constant(ngrid, ngrid, alphax_val); // row,col,value
  constexpr double alphay_val = 1e-10;
  Eigen::MatrixXd alphay = Eigen::MatrixXd::Constant(ngrid, ngrid, alphay_val); // row,col,value
  grid64.setAlpha(alphax, alphay);

  // create a boundary with constant values
  Boundary bc64 = Boundary(grid64);
  bc64.setBoundarySideClosed(BC_SIDE_LEFT);
  bc64.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc64.setBoundarySideClosed(BC_SIDE_TOP);
  bc64.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // set up a simulation environment
  Simulation Sim64 =
      Simulation(grid64, bc64, approach); // grid_64,boundary,simulation-approach

  //Sim64.setSolver(THOMAS_ALGORITHM_SOLVER);
 
  // set the timestep of the simulation
  Sim64.setTimestep(dt); // timestep

  // set the number of iterations
  Sim64.setIterations(2);

  // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  Sim64.setOutputCSV(CSV_OUTPUT_ON);

  // set output to the console to 'ON'
  Sim64.setOutputConsole(CONSOLE_OUTPUT_OFF);

  // // **** RUN SIM64 ****
  auto begin64 = std::chrono::high_resolution_clock::now();

  Sim64.run();

  auto end64 = std::chrono::high_resolution_clock::now();
  auto ms64 = std::chrono::duration_cast<std::chrono::milliseconds>(end64 - begin64);

  const double sum_after64 = grid64.getConcentrations().sum();

  std::cout << "Sum of init field:      " << std::setprecision(15) << sum_init64 
            << "\nSum after 2 iterations: " << sum_after64 
            << "\nMilliseconds: " << ms64.count() << std::endl << std::endl;
  return 0;
  
}

int DoFloat(int ngrid, APPROACH approach) {
  
  constexpr double dt = 10;
  
  // create a grid
  std::cout << "FLOAT grid: " << ngrid << std::endl;

  Grid32 grid32(ngrid, ngrid);

  // (optional) set the domain, e.g.:
  grid32.setDomain(0.1, 0.1);

  Eigen::MatrixXf initConc32 = Eigen::MatrixXf::Constant(ngrid, ngrid, 0);
  initConc32(50, 50) = 1E-6;
  grid32.setConcentrations(initConc32);
  
  const float sum_init32 = initConc32.sum();

  constexpr float alphax_val = 5e-10;
  Eigen::MatrixXf alphax = Eigen::MatrixXf::Constant(ngrid, ngrid, alphax_val); // row,col,value
  constexpr float alphay_val = 1e-10;
  Eigen::MatrixXf alphay = Eigen::MatrixXf::Constant(ngrid, ngrid, alphay_val); // row,col,value
  grid32.setAlpha(alphax, alphay);

  // create a boundary with constant values
  Boundary bc32 = Boundary(grid32);
  bc32.setBoundarySideClosed(BC_SIDE_LEFT);
  bc32.setBoundarySideClosed(BC_SIDE_RIGHT);
  bc32.setBoundarySideClosed(BC_SIDE_TOP);
  bc32.setBoundarySideClosed(BC_SIDE_BOTTOM);

  // set up a simulation environment
  Simulation Sim32 =
      Simulation(grid32, bc32, approach); // grid_32,boundary,simulation-approach

  // Sim32.setSolver(THOMAS_ALGORITHM_SOLVER);
 
  // set the timestep of the simulation
  Sim32.setTimestep(dt); // timestep

  // set the number of iterations
  Sim32.setIterations(2);

  // set kind of output [CSV_OUTPUT_OFF (default), CSV_OUTPUT_ON,
  // CSV_OUTPUT_VERBOSE]
  Sim32.setOutputCSV(CSV_OUTPUT_ON);

  // set output to the console to 'ON'
  Sim32.setOutputConsole(CONSOLE_OUTPUT_OFF);

  // // **** RUN SIM32 ****
  auto begin32 = std::chrono::high_resolution_clock::now();

  Sim32.run();

  auto end32 = std::chrono::high_resolution_clock::now();
  auto ms32 = std::chrono::duration_cast<std::chrono::milliseconds>(end32 - begin32);

  const float sum_after32 = grid32.getConcentrations().sum();

  std::cout << "Sum of init field:      " << std::setprecision(15) << sum_init32 
            << "\nSum after 2 iteration: " << sum_after32 
            << "\nMilliseconds: " << ms32.count() << std::endl << std::endl;
  return 0;
  
}


int main(int argc, char *argv[]) {

  int n[] = {101, 201, 501, 1001, 2001};

  for (int i = 0; i < std::size(n); i++) {
    DoFloat(n[i], BTCS_APPROACH);
    DoDble(n[i],  BTCS_APPROACH);
    DoFloat(n[i], FTCS_APPROACH);
    DoDble(n[i],  FTCS_APPROACH);
  }

  return 0;
}
