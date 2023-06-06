#include "tug/BoundaryCondition.hpp"
#include <tug/Diffusion.hpp>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace tug::diffusion;
using namespace tug::bc;

int main(int argc, char *argv[]) {

    // define problem dimensionality
    // set grid sizes for each dimension
    int dim = 1;
    int n = 20;

    vector<double> alpha(n, 1e-1);
    // double alpha = 1e-1;
    vector<double> field(n, 1e-6);
    for (int i = 1; i<20; i++) {
        field[i] = 0;
    }
    // double field = 1e-6;

    
    TugGrid grid_param; // why is grid_param defined separately?
    // grid_param.grid_cells[0] = 20;
    // grid_param.grid_cells[1] = 0;
    // grid_param.grid_cells[2] = 0;

    // grid_param.domain_size[0] = 20;
    // grid_param.domain_size[1] = 0;
    // grid_param.domain_size[2] = 0;

    

    TugInput input_param;
    input_param.setTimestep(1.);
    //input_param.grid = grid_param; 
    input_param.setGridCellN(n);
    input_param.setDomainSize(n); // what is domain?????
    BoundaryCondition bc(n);
    input_param.setBoundaryCondition(bc);

    BoundaryCondition bc2 = input_param.getBoundaryCondition();
    auto [bc_left, bc_right] = bc2.row_boundary(0);
    cout << "left: " << unsigned(bc_left.type) << endl;
    cout << "right: " << unsigned(bc_right.type) << endl;

    ofstream myfile;
    myfile.open("output.csv"); 
    if (!myfile) {
        exit(1);
    }

    for (int t = 0; t < 10000; t++) {
        double result = BTCS_1D(input_param, &field[0], &alpha[0]);
        //myfile << result;
        //myfile << '\n';
        //myfile << "Vector contents: ";
        for (int i = 0; i < field.size(); i++) {
            myfile << field[i];
            if (i < field.size()-1) {
                myfile << ", ";
            }
        }
        myfile << std::endl;
    }
    myfile.close();
}