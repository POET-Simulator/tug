#include <tug/Boundary.hpp>
#include <stdexcept>

using namespace std;

Boundary::Boundary(Grid &grid, BC_TYPE type) : grid(grid) {
    this->type = type;

    if (type == BC_TYPE_CONSTANT) {
        if (grid.getDim() == 1) {
            this->left = VectorXd::Constant(1, 1, 1);
            this->right = VectorXd::Constant(1, 1, 1);
        } else if (grid.getDim() == 2) {
            this->left = VectorXd::Constant(1, 1, 1);
            this->right = VectorXd::Constant(1, 1, 1);
        } else {
            throw invalid_argument("Dimension must be 1 or 2!");
        }
    }
}