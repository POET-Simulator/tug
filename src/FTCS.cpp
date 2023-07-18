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
                + timestep / (deltaRow*deltaRow) * (
                    calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
                    * grid.getConcentrations()(row+1,col)
                    - (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                    + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col)))
                    * grid.getConcentrations()(row,col)
                    + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
                    * grid.getConcentrations()(row-1,col)
                )
                - timestep / (deltaCol*deltaCol) * (
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
    int col = 0;
    for (int row = 1; row < rowMax-1; row++) {
        concentrations_t1(row, col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
                * grid.getConcentrations()(row,col+1)
                - (calc_alpha_intercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
                + 2 * grid.getAlphaX()(row,col)) * grid.getConcentrations()(row,col) 
                + 2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_LEFT)(row, 1))
            + timestep / (deltaRow*deltaRow)
                * (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
                * grid.getConcentrations()(row+1,col)
                - (calc_alpha_intercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col)))
                * grid.getConcentrations()(row,col)
                + calc_alpha_intercell(grid.getAlphaY()(row-1,col), grid.getConcentrations()(row,col))
                * grid.getConcentrations()(row-1,col));
    }

    // right without corners / looping over columns
    col = colMax-1;
    for (int row = 1; row < rowMax-1; row++) {
        concentrations_t1(row,col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_RIGHT)(row, 1)
                - (calc_alpha_intercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
                + 2 * grid.getAlphaX()(row,col)) + 2 * grid.getAlphaX()(row,col) 
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
    for(int col=1; col<colMax-1;col++){
        int row = 0;
        concentrations_t1(row, col) = grid.getConcentrations()(row, col)
            + timestep/(grid.getDeltaRow()*grid.getDeltaRow()) * (calc_alpha_intercell(grid.getAlphaY()(1, col), grid.getAlphaY()(0, col)) * grid.getConcentrations()(1,col)
                - (calc_alpha_intercell(grid.getAlphaY()(1, col), grid.getAlphaY()(0, col)) + 2 * grid.getAlphaY()(0, col)) * grid.getConcentrations()(0, col)
                + 2 * grid.getAlphaY()(0, col) * bc.getBoundaryConditionValue(BC_SIDE_TOP)(1, col))
            + timestep/(grid.getDeltaCol()*grid.getDeltaCol()) * (calc_alpha_intercell(grid.getAlphaX()(0, col+1), grid.getAlphaX()(0, col)) * grid.getConcentrations()(0, col+1)
                - (calc_alpha_intercell(grid.getAlphaX()(0, col+1), grid.getAlphaX()(0, col)) + calc_alpha_intercell(grid.getAlphaX()(0, col-1), grid.getAlphaX()(0, col))) * grid.getConcentrations()(0, col)
                + calc_alpha_intercell(grid.getAlphaX()(0, col-1), grid.getAlphaX()(0, col)) * grid.getConcentrations()(0, col-1));
    }

    // bottom without corners / looping over cols
    int row = rowMax-1;
    for(int col=1; row<colMax-1;col++){
        concentrations_t1(row, col) = grid.getConcentrations()(row, col)
            + timestep/(grid.getDeltaRow()*grid.getDeltaRow()) * (2 * grid.getAlphaY()(row, col) * bc.getBoundaryConditionValue(BC_SIDE_BOTTOM)(1, col)
                - (calc_alpha_intercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)) + 2 * grid.getAlphaY()(row, col)) * grid.getConcentrations()(row, col)
                + calc_alpha_intercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)))
            + timestep/(grid.getDeltaCol()*grid.getDeltaCol()) * (calc_alpha_intercell(grid.getAlphaX()(row, col+1), grid.getAlphaX()(row, col)) * grid.getConcentrations()(row, col+1)
                - (calc_alpha_intercell(grid.getAlphaX()(row, col+1), grid.getAlphaX()(row, col)) + calc_alpha_intercell(grid.getAlphaX()(row, col-1), grid.getAlphaX()(row, col))) * grid.getConcentrations()(row, col)
                + calc_alpha_intercell(grid.getAlphaX()(row, col-1), grid.getAlphaX()(row, col-1)) * grid.getConcentrations()(row, col-1));
    }


    concentrations_t1(0,0) = 0;
    concentrations_t1(rowMax-1,0) = 0;
    concentrations_t1(0,colMax-1) = 0;
    concentrations_t1(rowMax-1,colMax-1) = 0;

    grid.setConcentrations(concentrations_t1);
}

auto FTCS_closed(Grid &grid, Boundary &bc, double timestep) {
    return 0;
}

void FTCS(Grid &grid, Boundary &bc, double timestep) {
    if (bc.getBoundaryConditionType() == BC_TYPE_CONSTANT) {
        FTCS_constant(grid, bc, timestep);
    } else if (bc.getBoundaryConditionType() == BC_TYPE_CLOSED) {
        FTCS_closed(grid, bc, timestep);
    }
}
