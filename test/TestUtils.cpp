#include <ios>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace Eigen;

MatrixXd CSV2Eigen(string file2Convert) {

    vector<double> matrixEntries;

    ifstream matrixDataFile(file2Convert);
    if (matrixDataFile.fail()) {
        throw invalid_argument("File probably non-existent!");
    }

    string matrixRowString;
    string matrixEntry;
    int matrixRowNumber = 0;

    while(getline(matrixDataFile, matrixRowString)){
        stringstream matrixRowStringStream(matrixRowString);
        while(getline(matrixRowStringStream, matrixEntry, ' ')){
            matrixEntries.push_back(stod(matrixEntry));
        }
        if (matrixRowString.length() > 1) {
            matrixRowNumber++;
        }
    }

    return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}

bool checkSimilarity(MatrixXd a, MatrixXd b, double precision=1e-5) {
    return a.isApprox(b, precision);
}

bool checkSimilarityV2(MatrixXd a, MatrixXd b, double maxDiff) {

    MatrixXd diff = a - b;
    double maxCoeff = diff.maxCoeff();
    return abs(maxCoeff) < maxDiff;
}