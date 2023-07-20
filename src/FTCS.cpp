#include <cstddef>
#include <tug/Boundary.hpp>
#include <iostream>

using namespace std;

double calc_alpha_intercell(double alpha1, double alpha2, bool useHarmonic = false) {
    if (useHarmonic) {
        return 2 / ((1/alpha1) + (1/alpha2));
    } else {
        return 0.5 * (alpha1 + alpha2);
    }
}

// IN PROGRESS
double verticalChange(Grid grid, int row, int col) {
    
    double result =    
                calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
                * grid.getConcentrations()(row+1,col)
                - (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col)))
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
                * grid.getConcentrations()(row-1,col);

    return result;
}

// IN PROGRESS
double horizontalChange(Grid grid, int row, int col) {
    
    double result = 
                calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
                * grid.getConcentrations()(row,col+1)
                - (calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col))
                + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col)))
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                * grid.getConcentrations()(row,col-1);

    return result;
}

// IN PROGRESS
MatrixXd FTCS_1D_constant(Grid grid, Boundary bc, double timestep) {
    int colMax = grid.getCol();
    double deltaCol = grid.getDeltaCol();

    MatrixXd concentrations_t1 = MatrixXd::Constant(1, colMax, 0);

    int row = 0;
    // for (int col = 1; col < colMax-1; col++) {
    //     concentrations_t1 = grid.getConcentrations()(row,col)
    //         + horizontal_term();
    // }
}

// IN PROGRESS
MatrixXd FTCS_2D_constant(Grid grid, Boundary bc, double timestep) {
    int rowMax = grid.getRow();
    int colMax = grid.getCol();
    double deltaRow = grid.getDeltaRow();
    double deltaCol = grid.getDeltaCol();

}

MatrixXd FTCS_constant(Grid grid, Boundary bc, double timestep) {
    int rowMax = grid.getRow();
    int colMax = grid.getCol();
    double deltaRow = grid.getDeltaRow();
    double deltaCol = grid.getDeltaCol();

    // Matrix with concentrations at time t+1
    // TODO profiler / only use 2 matrices
    MatrixXd concentrations_t1 = MatrixXd::Constant(rowMax, colMax, 0);

    // inner cells
    // (should have 7 calls to current concentration)
    for (int row = 1; row < rowMax-1; row++) {
        for (int col = 1; col < colMax-1; col++) {
            concentrations_t1(row, col) = grid.getConcentrations()(row, col) 
                + timestep / (deltaRow*deltaRow) 
                    * (
                    calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
                    * grid.getConcentrations()(row+1,col)
                    - (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                    + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col)))
                    * grid.getConcentrations()(row,col)
                    + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
                    * grid.getConcentrations()(row-1,col)
                    )
                + timestep / (deltaCol*deltaCol) 
                    * (
                    calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
                    * grid.getConcentrations()(row,col+1)
                    - (calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col))
                    + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col)))
                    * grid.getConcentrations()(row,col)
                    + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                    * grid.getConcentrations()(row,col-1)
                    )
                ;
        }
    }

    // boundary conditions
    // left without corners / looping over rows
    // (should have 6 calls to current concentration)
    int col = 0;
    for (int row = 1; row < rowMax-1; row++) {
        concentrations_t1(row, col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
                * grid.getConcentrations()(row,col+1)
                - (calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
                + 2 * grid.getAlphaX()(row,col)) * grid.getConcentrations()(row,col) 
                + 2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_LEFT)(row))
            + timestep / (deltaRow*deltaRow)
                * (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
                * grid.getConcentrations()(row+1,col)
                - (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col)))
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
                * grid.getConcentrations()(row-1,col));
    }

    // right without corners / looping over columns
    // (should have 6 calls to current concentration)
    col = colMax-1;
    for (int row = 1; row < rowMax-1; row++) {
        concentrations_t1(row,col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_RIGHT)(row)
                - (calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                + 2 * grid.getAlphaX()(row,col)) 
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                * grid.getConcentrations()(row,col-1))
            + timestep / (deltaRow*deltaRow)
                * (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                * grid.getConcentrations()(row+1,col)
                - (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col)))
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
                * grid.getConcentrations()(row-1,col));
    }


    // top without corners / looping over cols
    // (should have 6 calls to current concentration)
    int row = 0;
    for (int col=1; col<colMax-1;col++){
        concentrations_t1(row, col) = grid.getConcentrations()(row, col)
            + timestep/(grid.getDeltaRow()*grid.getDeltaRow()) * (calc_alpha_intercell(grid.getAlphaY()(1, col), grid.getAlphaY()(0, col)) * grid.getConcentrations()(1,col)
                - (calc_alpha_intercell(grid.getAlphaY()(1, col), grid.getAlphaY()(0, col)) + 2 * grid.getAlphaY()(0, col)) * grid.getConcentrations()(0, col)
                + 2 * grid.getAlphaY()(0, col) * bc.getBoundaryConditionValue(BC_SIDE_TOP)(col))
            + timestep/(grid.getDeltaCol()*grid.getDeltaCol()) * (calc_alpha_intercell(grid.getAlphaX()(0, col+1), grid.getAlphaX()(0, col)) * grid.getConcentrations()(0, col+1)
                - (calc_alpha_intercell(grid.getAlphaX()(0, col+1), grid.getAlphaX()(0, col)) + calc_alpha_intercell(grid.getAlphaX()(0, col-1), grid.getAlphaX()(0, col))) * grid.getConcentrations()(0, col)
                + calc_alpha_intercell(grid.getAlphaX()(0, col-1), grid.getAlphaX()(0, col)) * grid.getConcentrations()(0, col-1));
    }

    // bottom without corners / looping over cols
    // (should have 6 calls to current concentration)
    row = rowMax-1;
    for(int col=1; col<colMax-1;col++){
        concentrations_t1(row, col) = grid.getConcentrations()(row, col)
            + timestep/(grid.getDeltaRow()*grid.getDeltaRow()) * (2 * grid.getAlphaY()(row, col) * bc.getBoundaryConditionValue(BC_SIDE_BOTTOM)(col)
                - (calc_alpha_intercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)) + 2 * grid.getAlphaY()(row, col)) * grid.getConcentrations()(row, col)
                + calc_alpha_intercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)) * grid.getConcentrations()(row-1,col))
            + timestep/(grid.getDeltaCol()*grid.getDeltaCol()) * (calc_alpha_intercell(grid.getAlphaX()(row, col+1), grid.getAlphaX()(row, col)) * grid.getConcentrations()(row, col+1)
                - (calc_alpha_intercell(grid.getAlphaX()(row, col+1), grid.getAlphaX()(row, col)) + calc_alpha_intercell(grid.getAlphaX()(row, col-1), grid.getAlphaX()(row, col))) * grid.getConcentrations()(row, col)
                + calc_alpha_intercell(grid.getAlphaX()(row, col-1), grid.getAlphaX()(row, col)) * grid.getConcentrations()(row, col-1));
    }

    // corner top left
    // (should have 5 calls to current concentration)
    row = 0;
    col = 0;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col)
        + timestep/(deltaCol*deltaCol)
            * (
                calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) * grid.getConcentrations()(row,col+1)
                - ( 
                    calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col))
                    + 2 * grid.getAlphaX()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + 2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_LEFT)(row)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) * grid.getConcentrations()(row+1,col)
                - (
                    calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                    + 2 * grid.getAlphaY()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + 2 * grid.getAlphaY()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_TOP)(col)
            )
        ;

    // corner top right
    // (should have 5 calls to current concentration)
    row = 0;
    col = colMax-1;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col) 
        + timestep/(deltaCol*deltaCol)
            * (
                2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_RIGHT)(row)
                - (
                    calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                    + 2 * grid.getAlphaX()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                * grid.getConcentrations()(row,col-1)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) * grid.getConcentrations()(row+1,col)
                - (
                    calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                    + 2 * grid.getAlphaY()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + 2 * grid.getAlphaY()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_TOP)(col)
            )
        ;

    // corner bottom left
    // (should have 5 calls to current concentration)
    row = rowMax-1;
    col = 0;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col)
        + timestep/(deltaCol*deltaCol)
            * (
                calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) * grid.getConcentrations()(row,col+1)
                - ( 
                    calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col))
                    + 2 * grid.getAlphaX()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + 2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_LEFT)(row)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                2 * grid.getAlphaY()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_BOTTOM)(col)
                - (
                    calc_alpha_intercell(grid.getAlphaY()(row,col), grid.getAlphaY()(row-1,col))
                    + 2 * grid.getAlphaY()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaY()(row,col), grid.getAlphaY()(row-1,col))
                * grid.getConcentrations()(row-1,col)
            )
        ;

    // corner bottom right
    // (should have 5 calls to current concentration)
    row = rowMax-1;
    col = colMax-1;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col) 
        + timestep/(deltaCol*deltaCol)
            * (
                2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_RIGHT)(row)
                - (
                    calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                    + 2 * grid.getAlphaX()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                * grid.getConcentrations()(row,col-1)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                2 * grid.getAlphaY()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_BOTTOM)(col)
                - (
                    calc_alpha_intercell(grid.getAlphaY()(row,col), grid.getAlphaY()(row-1,col))
                    + 2 * grid.getAlphaY()(row,col)
                )
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaY()(row,col), grid.getAlphaY()(row-1,col))
                * grid.getConcentrations()(row-1,col)
            )
        ;


    return concentrations_t1;
}

// TODO
MatrixXd FTCS_closed(Grid grid, Boundary bc, double timestep) {
    return MatrixXd();
}

MatrixXd FTCS(Grid grid, Boundary bc, double timestep) {
    switch (bc.getBoundaryConditionType()) {
        case BC_TYPE_CONSTANT:
            return FTCS_constant(grid, bc, timestep);
        case BC_TYPE_CLOSED:
            return FTCS_closed(grid, bc, timestep);
        default:
            // TODO handle
            return MatrixXd();
    }
}
