#include "tug/BoundaryCondition.hpp"
#include <tug/Diffusion.hpp>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace tug::diffusion;
using namespace tug::bc;

int main(int argc, char *argv[]) {

    int dim = 2;
    int n = 20;
    int m = 20;

    vector<double> alpha(n * m, 1);
    vector<double> field(n * m, 0);
    field[n * 19] = 2000;
    field[n * 19 + 19] = 2000;

    // for (int i = 1; i<20; i++) {
    //     for (int j = 0; j<20; j++ ) {
    //         field[i] = 0;
    //     }
    // }

    // print field
    cout << "Initial field:" << endl;
    for (int i = 0; i<n; i++) {
        for (int j = 0; j<m; j++) {
            cout << field[n * i + j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    TugInput input_param;
    input_param.setTimestep(0.1);
    input_param.setGridCellN(n, m);
    input_param.setDomainSize(n, m);
    BoundaryCondition bc(n, m);
    boundary_condition bc_constant;
    bc_constant.type = BC_TYPE_CLOSED;
    bc_constant.value = 0;
    bc.setSide(BC_SIDE_LEFT, bc_constant);
    bc.setSide(BC_SIDE_RIGHT, bc_constant);
    bc.setSide(BC_SIDE_BOTTOM, bc_constant);
    // bc_constant.value = 2000;
    bc.setSide(BC_SIDE_TOP, bc_constant);
    input_param.setBoundaryCondition(bc);

    // int iterations = 1000;
    // for (int t = 0; t < iterations; t++) {
    //     double result = ADI_2D(input_param, &field[0], &alpha[0]);

    //     if (t % 100 == 0) {
    //         cout << "Iteration " << t << ":" << endl;
    //         for (int i = 0; i<n; i++) {
    //             for (int j = 0; j<m; j++) {
    //                 cout << field[n * i + j] << " ";
    //             }
    //             cout << endl;
    //         }
    //         cout << endl;
    //     }
    // }

    ofstream myfile;
    myfile.open("output.csv"); 
    if (!myfile) {
        exit(1);
    }

    int iterations = 1000;
    for (int t = 0; t < iterations; t++) {
        double result = ADI_2D(input_param, &field[0], &alpha[0]);
    
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                myfile << field[i * m + j];
                if (j < m-1) {
                    myfile << ", ";
                }
            }
            myfile << "\n";
        }
        myfile << std::endl << std::endl;
    }
    myfile.close();
}