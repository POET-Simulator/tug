#include "tug/BoundaryCondition.hpp"
#include <tug/Diffusion.hpp>

#include <iostream>

using namespace std;
using namespace tug::diffusion;
using namespace tug::bc;

int main(int argc, char *argv[]) {

    int dim = 2;
    int n = 5;
    int m = 5;

    vector<double> alpha(n * m, 1e-1);
    vector<double> field(n * m, 0);
    field[0] = 1e-6;
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
    input_param.setTimestep(1.);
    input_param.setGridCellN(n, m);
    input_param.setDomainSize(n, m);
    BoundaryCondition bc(n, m);
    input_param.setBoundaryCondition(bc);

    int iterations = 1000;
    for (int t = 0; t < iterations; t++) {
        double result = ADI_2D(input_param, &field[0], &alpha[0]);

        if (t % 100 == 0) {
            cout << "Iteration " << t << ":" << endl;
            for (int i = 0; i<n; i++) {
                for (int j = 0; j<m; j++) {
                    cout << field[n * i + j] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }
}