#include <algorithm>
#include <tug/BoundaryCondition.hpp>
#include <tug/Diffusion.hpp>
#include <vector>
#include <iostream>

using namespace std;
using namespace tug::bc;
using namespace tug::diffusion;

int main(int argc, char *argv[]) {
    TugInput input;

    BoundaryCondition example(10);

    uint8_t side = BC_SIDE_LEFT;
    boundary_condition bc;
    bc.type = BC_TYPE_CONSTANT;
    bc.value = 1;
    input.setBoundaryCondition(example);
    BoundaryCondition returnedBC = input.getBoundaryCondition();

    // example.setSide(side, bc);
    vector<boundary_condition> result_left = example.getSide(side);
    // vector<boundary_condition> result_top = example.getSide(BC_SIDE_TOP);

    for (auto i : result_left) {
        cout << i.value << " ";
    }
    cout << endl;
    // for (auto i : result_top) {
    //     cout << i.value << " ";
    // }
    // cout << endl;
}