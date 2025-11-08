#ifndef _IO_H
#define _IO_H

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <vector>

#include <Eigen/Eigen>

static std::string header;

constexpr int CSV_OUT_PRECISION = 18;

template <typename T> inline T convert(const std::string &str) {
  std::istringstream ss(str);
  T value;
  ss >> value;
  return value;
}

// read a csv file into Eigen::Matrix
template <class T>
inline std::vector<std::vector<T>> read_conc_csv(const std::string &filename,
                                                 std::size_t cols,
                                                 std::size_t nrows) {
  std::ifstream in(filename);
  if (!in.is_open())
    throw std::runtime_error("Could not open file");

  std::string line, field;

  // store header globally
  std::getline(in, header);

  std::size_t val_count = 0;
  std::stringstream head_ss(header);

  while (std::getline(head_ss, field, ',')) {
    val_count++;
  }

  std::vector<std::vector<T>> data(val_count);

  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::size_t colIdx = 0;
    while (std::getline(ss, field, ',')) {
      data[colIdx++].push_back(convert<T>(field));
    }
  }

  // std::vector<Eigen::MatrixX<T>> result(data.size());

  // for (int i = 0; i < data.size(); i++) {
  //   result[i] = Eigen::Map<
  //       Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
  //       data[i].data(), nrows, cols);
  // }
  return data;
}

template <class T>
inline std::vector<T> read_alpha_csv(const std::string &filename) {
  std::ifstream in(filename);
  if (!in.is_open())
    throw std::runtime_error("Could not open file");

  std::string line, field;

  std::vector<T> data;

  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::size_t colIdx = 0;
    while (std::getline(ss, field, ',')) {
      const T value = convert<T>(field);
      data.push_back(value);
    }
  }

  return data;
}

// write an Eigen::Matrix to a csv file
template <class T>
inline void write_conc_csv(const std::string &filename,
                           const std::vector<std::vector<T>> &mat_strg) {
  std::ofstream out(filename);
  if (!out.is_open())
    throw std::runtime_error("Could not open file");

  out << header << std::endl;

  std::size_t rows = mat_strg[0].size();

  for (std::size_t i_row = 0; i_row < rows; i_row++) {
    for (std::size_t i_col = 0; i_col < mat_strg.size(); i_col++) {
      out << std::scientific << std::setprecision(CSV_OUT_PRECISION)
          << mat_strg[i_col][i_row];
      if (i_col < mat_strg.size() - 1)
        out << ",";
    }
    out << std::endl;
  }

  out.close();
}

#endif // _IO_H