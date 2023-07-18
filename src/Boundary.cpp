#include <iostream>
#include <tug/Boundary.hpp>
#include <stdexcept>

using namespace std;

Boundary::Boundary(Grid &grid, BC_TYPE type) : grid(grid) {
    //probably to DEBUG assignment grid
    this->type = type;

    if (type == BC_TYPE_CONSTANT) {
        if (grid.getDim() == 1) {
            this->left = VectorXd::Constant(1, 1);
            this->right = VectorXd::Constant(1, 1);
        } else if (grid.getDim() == 2) {
            this->left = VectorXd::Constant(grid.getRow(), 1);
            this->right = VectorXd::Constant(grid.getRow(), 1);
            this->top = VectorXd::Constant(grid.getCol(), 1);
            this->bottom = VectorXd::Constant(grid.getCol(), 1);

        } else {
            throw invalid_argument("Dimension must be 1 or 2!");
        }
    } 
}

BC_TYPE Boundary::getBoundaryConditionType() {
    return this->type; 
}

void Boundary::setBoundaryConditionValue(BC_SIDE side, VectorXd &values) {
    if (type != BC_TYPE_CONSTANT) {
        // TODO check if correct way for handling warning
        cerr << "Values will not be used, wrong BC_TYPE!";
    }

    switch (side) {
        case BC_SIDE_LEFT:
            this->left = values;
            break;
        case BC_SIDE_RIGHT:
            this->right = values;
            break;
        case BC_SIDE_TOP:
            this->top = values;
            break;
        case BC_SIDE_BOTTOM:
            this->bottom = values;
            break;
        default:
            throw invalid_argument("Invalid side given!");
    }
}

VectorXd Boundary::getBoundaryConditionValue(BC_SIDE side) {
    switch (side) {
        case BC_SIDE_LEFT:
            return this->left;
        case BC_SIDE_RIGHT:
            return this->right;
        case BC_SIDE_TOP:
            return this->top;
        case BC_SIDE_BOTTOM:
            return this->bottom;
        default:
            throw invalid_argument("Invalid side given!");
    }
}