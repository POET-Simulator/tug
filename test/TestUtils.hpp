#include "tug/Core/Matrix.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <gtest/gtest.h>

#define TUG_TEST(x) TEST(Tug, x)

inline tug::RowMajMat<double> CSV2Eigen(std::string file2Convert) {

  std::vector<double> matrixEntries;

  std::ifstream matrixDataFile(file2Convert);
  if (matrixDataFile.fail()) {
    throw std::invalid_argument("File probably non-existent!");
  }

  std::string matrixRowString;
  std::string matrixEntry;
  int matrixRowNumber = 0;

  while (getline(matrixDataFile, matrixRowString)) {
    std::stringstream matrixRowStringStream(matrixRowString);
    while (getline(matrixRowStringStream, matrixEntry, ' ')) {
      matrixEntries.push_back(stod(matrixEntry));
    }
    if (matrixRowString.length() > 1) {
      matrixRowNumber++;
    }
  }

  return tug::RowMajMatMap<double>(matrixEntries.data(), matrixRowNumber,
                                   matrixEntries.size() / matrixRowNumber);
}

inline bool checkSimilarity(tug::RowMajMat<double> &a,
                            tug::RowMajMatMap<double> &b,
                            double precision = 1e-5) {
  return a.isApprox(b, precision);
}

inline bool checkSimilarityV2(tug::RowMajMat<double> &a,
                              tug::RowMajMatMap<double> &b, double maxDiff) {

  tug::RowMajMat<double> diff = a - b;
  double maxCoeff = diff.maxCoeff();
  return abs(maxCoeff) < maxDiff;
}
