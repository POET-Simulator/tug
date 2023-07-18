#include <tug/Boundary.hpp>


auto calc_alpha_intercell(double alpha1, double alpha2, bool useHarmonic = false) {
    if (useHarmonic) {
        return 2 / ((1/alpha1) + (1/alpha2));
    } else {
        return 0.5 * (alpha1 + alpha2);
    }
}

auto FTCS_constant(Grid &grid, Boundary &bc, double timestep) {
    int rowMax = grid.getRow();
    int colMax = grid.getCol();
    double deltaRow = grid.getDeltaRow();
    double deltaCol = grid.getDeltaCol();

    // Matrix with concentrations at time t+1
    // TODO profiler / only use 2 matrices
    Matrix2d concentrations_t1 = Matrix2d(rowMax, colMax);

    // inner cells
    for (int row = 1; row < rowMax-1; row++) {
        for (int col = 1; col < colMax-1; col++) {
            concentrations_t1(row, col) = grid.getConcentrations()(row, col) 
                + timestep / deltaRow*deltaRow * (
                    calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
                    * grid.getConcentrations()(row+1,col)
                    - (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                    + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col)))
                    * grid.getConcentrations()(row,col)
                    + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
                    * grid.getConcentrations()(row-1,col)
                )
                - timestep / deltaCol*deltaCol * (
                    calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
                    * grid.getConcentrations()(row,col+1)
                    - (calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col))
                    + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col)))
                    * grid.getConcentrations()(row,col)
                    + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                    * grid.getConcentrations()(row,col-1)
                );
        }
    }

    // boundary conditions
    // left without corners / looping over rows

    return 0;
}

auto FTCS_closed(Grid &grid, Boundary &bc, double timestep) {
    return 0;
}

auto FTCS(Grid &grid, Boundary &bc, double timestep) {
    if (bc.getBoundaryConditionType() == BC_TYPE_CONSTANT) {
        int test = FTCS_constant(grid, bc, timestep);
    } else if (bc.getBoundaryConditionType() == BC_TYPE_CLOSED) {
        FTCS_closed(grid, bc, timestep);
    }
}
