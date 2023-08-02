#include <iostream>
#include <Eigen/Core>
#include <fstream>

using namespace std;
using namespace Eigen;

inline MatrixXd CSV2Eigen(string file2Convert){

    vector<double> matrixEntries;

    ifstream matrixDataFile(file2Convert);

    string matrixRowString;
    string matrixEntry;
    int matrixRowNumber = 0;

    while(getline(matrixDataFile, matrixRowString)){
        stringstream matrixRowStringStream(matrixRowString);
        while(getline(matrixRowStringStream, matrixEntry, ' ')){
            matrixEntries.push_back(stod(matrixEntry));
        }
        matrixRowNumber++;
    }

    return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);

}