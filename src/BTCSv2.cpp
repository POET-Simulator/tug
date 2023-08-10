#include "TugUtils.hpp"
#include <tug/Boundary.hpp>

using namespace Eigen;


// calculates arithmetic or harmonic mean of alpha between two cells
static double calcAlphaIntercell(double &alpha1, double &alpha2, bool useHarmonic = true) {
    if (useHarmonic) {
        return double(2) / ((double(1)/alpha1) + (double(1)/alpha2));
    } else {
        return 0.5 * (alpha1 + alpha2);
    }
}


static MatrixXd createCoeffMatrix(Grid &grid, int rowIndex, double sx) {

    // square matrix of column^2 dimension for the coefficients
    int dim = grid.getCol();
    SparseMatrix<double, RowMajor> cm(dim, dim);

    // top left
    cm.coeffRef(0,0) = 1 + sx * (calcAlphaIntercell(grid.getAlphaX()(rowIndex,0), grid.getAlphaX()(rowIndex,1)));


}


// BTCS solution for 1D grid
static void BTCS_1D(Grid &grid, Boundary &bc, double &timestep) {

}


// BTCS solution for 2D grid
static void BTCS_2D(Grid &grid, Boundary &bc, double &timestep) {

}


// entry point; differentiate between 1D and 2D grid
static void BTCS(Grid &grid, Boundary &bc, double &timestep) {
    if (grid.getDim() == 1) {
        BTCS_1D(grid, bc, timestep);
    } else if (grid.getDim() == 2) {
        BTCS_2D(grid, bc, timestep);
    } else {
        throw_invalid_argument("Error: Only 1- and 2-dimensional grids are defined!");
    }
}