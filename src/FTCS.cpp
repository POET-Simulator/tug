#include <cstddef>
#include <tug/Boundary.hpp>
#include <iostream>

using namespace std;

double calcAlphaIntercell(double alpha1, double alpha2, bool useHarmonic = false) {
    if (useHarmonic) {
        return 2 / ((1/alpha1) + (1/alpha2));
    } else {
        return 0.5 * (alpha1 + alpha2);
    }
}


double calcHorizontalChange(Grid grid, int row, int col) {

    double result = 
        calcAlphaIntercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
            * grid.getConcentrations()(row,col+1)
        - (
            calcAlphaIntercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col))
            + calcAlphaIntercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
            )
            * grid.getConcentrations()(row,col)
        + calcAlphaIntercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
            * grid.getConcentrations()(row,col-1);

    return result;
}


double calcVerticalChange(Grid grid, int row, int col) {
    
    double result =    
        calcAlphaIntercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col)) 
            * grid.getConcentrations()(row+1,col)
        - (
            calcAlphaIntercell(grid.getAlphaY()(row+1,col), grid.getAlphaY()(row,col))
            + calcAlphaIntercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
            )
            * grid.getConcentrations()(row,col)
        + calcAlphaIntercell(grid.getAlphaY()(row-1,col), grid.getAlphaY()(row,col))
            * grid.getConcentrations()(row-1,col);

    return result;
}


double calcHorizontalChangeLeftBoundary(Grid grid, Boundary bc, int row, int col) {

    double result = 
        calcAlphaIntercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
            * grid.getConcentrations()(row,col+1)
        - (
            calcAlphaIntercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
            + 2 * grid.getAlphaX()(row,col)
            ) 
            * grid.getConcentrations()(row,col) 
        + 2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_LEFT)(row);

    return result;
}


double calcHorizontalChangeRightBoundary(Grid grid, Boundary bc, int row, int col) {

    double result = 
        2 * grid.getAlphaX()(row,col) * bc.getBoundaryConditionValue(BC_SIDE_RIGHT)(row)
        - (
            calcAlphaIntercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
            + 2 * grid.getAlphaX()(row,col)
            ) 
            * grid.getConcentrations()(row,col)
        + calcAlphaIntercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
            * grid.getConcentrations()(row,col-1);

    return result;
}


double calcVerticalChangeTopBoundary(Grid grid, Boundary bc, int row, int col) {
    
    double result = 
        calcAlphaIntercell(grid.getAlphaY()(row+1, col), grid.getAlphaY()(row, col)) 
            * grid.getConcentrations()(row+1,col)
        - (
            calcAlphaIntercell(grid.getAlphaY()(row+1, col), grid.getAlphaY()(row, col)) 
            + 2 * grid.getAlphaY()(row, col)
            ) 
            * grid.getConcentrations()(row, col)
        + 2 * grid.getAlphaY()(row, col) * bc.getBoundaryConditionValue(BC_SIDE_TOP)(col);

    return result;
}


double calcVerticalChangeBottomBoundary(Grid grid, Boundary bc, int row, int col) {

    double result = 
        2 * grid.getAlphaY()(row, col) * bc.getBoundaryConditionValue(BC_SIDE_BOTTOM)(col)
        - (
            calcAlphaIntercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)) 
            + 2 * grid.getAlphaY()(row, col)
            ) 
            * grid.getConcentrations()(row, col)
        + calcAlphaIntercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)) 
            * grid.getConcentrations()(row-1,col);

    return result;
}


MatrixXd FTCS_1D(Grid grid, Boundary bc, double timestep) {
    int colMax = grid.getCol();
    double deltaCol = grid.getDeltaCol();

    MatrixXd concentrations_t1 = MatrixXd::Constant(1, colMax, 0);

    // only one row in 1D case
    int row = 0;

    // inner cells
    for (int col = 1; col < colMax-1; col++) {
        concentrations_t1(row,col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol)
                * (
                    calcHorizontalChange(grid, row, col)
                )
            ;
    }

    // left boundary
    int col = 0;
    concentrations_t1(row, col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChangeLeftBoundary(grid, bc, row, col)
                )
            ;


    // right boundary
    col = colMax-1;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChangeRightBoundary(grid, bc, row, col)
                )
            ;

    return concentrations_t1;
}


MatrixXd FTCS_2D(Grid grid, Boundary bc, double timestep) {
    int rowMax = grid.getRow();
    int colMax = grid.getCol();
    double deltaRow = grid.getDeltaRow();
    double deltaCol = grid.getDeltaCol();

    // Matrix with concentrations at time t+1
    // TODO profiler / only use 2 matrices
    MatrixXd concentrations_t1 = MatrixXd::Constant(rowMax, colMax, 0);

    // inner cells
    // these do not depend on the boundary condition type
    for (int row = 1; row < rowMax-1; row++) {
        for (int col = 1; col < colMax-1; col++) {
            concentrations_t1(row, col) = grid.getConcentrations()(row, col) 
                + timestep / (deltaRow*deltaRow) 
                    * (
                        calcVerticalChange(grid, row, col)
                    )
                + timestep / (deltaCol*deltaCol) 
                    * (
                        calcHorizontalChange(grid, row, col)
                    )
                ;
        }
    }

    // boundary conditions
    // left without corners / looping over rows
    int col = 0;
    for (int row = 1; row < rowMax-1; row++) {
        concentrations_t1(row, col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChangeLeftBoundary(grid, bc, row, col)
                )
            + timestep / (deltaRow*deltaRow)
                * (
                    calcVerticalChange(grid, row, col)
                )
            ;
    }

    // right without corners / looping over columns
    col = colMax-1;
    for (int row = 1; row < rowMax-1; row++) {
        concentrations_t1(row,col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChangeRightBoundary(grid, bc, row, col)
                )
            + timestep / (deltaRow*deltaRow)
                * (
                    calcVerticalChange(grid, row, col)
                )
            ;
    }


    // top without corners / looping over cols
    int row = 0;
    for (int col=1; col<colMax-1;col++){
        concentrations_t1(row, col) = grid.getConcentrations()(row, col)
            + timestep / (deltaRow*deltaRow) 
                * (
                    calcVerticalChangeTopBoundary(grid, bc, row, col)
                )
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChange(grid, row, col)
                )
            ;
    }

    // bottom without corners / looping over cols
    row = rowMax-1;
    for(int col=1; col<colMax-1;col++){
        concentrations_t1(row, col) = grid.getConcentrations()(row, col)
            + timestep / (deltaRow*deltaRow) 
                * (
                    calcVerticalChangeBottomBoundary(grid, bc, row, col)
                )
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChange(grid, row, col)
                )
            ;
    }

    // corner top left
    row = 0;
    col = 0;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col)
        + timestep/(deltaCol*deltaCol)
            * (
                calcHorizontalChangeLeftBoundary(grid, bc, row, col)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                calcVerticalChangeTopBoundary(grid, bc, row, col)
            )
        ;

    // corner top right
    row = 0;
    col = colMax-1;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col) 
        + timestep/(deltaCol*deltaCol)
            * (
                calcHorizontalChangeRightBoundary(grid, bc, row, col)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                calcVerticalChangeTopBoundary(grid, bc, row, col)
            )
        ;

    // corner bottom left
    row = rowMax-1;
    col = 0;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col)
        + timestep/(deltaCol*deltaCol)
            * (
                calcHorizontalChangeLeftBoundary(grid, bc, row, col)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                calcVerticalChangeBottomBoundary(grid, bc, row, col)
            )
        ;

    // corner bottom right
    row = rowMax-1;
    col = colMax-1;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col) 
        + timestep/(deltaCol*deltaCol)
            * (
                calcHorizontalChangeRightBoundary(grid, bc, row, col)
            )
        + timestep/(deltaRow*deltaRow)
            * (
                calcVerticalChangeBottomBoundary(grid, bc, row, col)
            )
        ;


    return concentrations_t1;
}


MatrixXd FTCS(Grid grid, Boundary bc, double timestep) {
    // inner cells 
    // TODO only the boundary cells are different in constant and closed case

    // if 1D:
    //      do inner cells 
    //      do left boundary according to bc type
    //      do right boundary according to bc type
    // if 2D:
    //      do inner cells
    //      do left boundaries according to bc type
    //      do right boundaries according to bc type
    //      ...

    if (grid.getDim() == 1) {
        return FTCS_1D(grid, bc, timestep);
    } else {
        return FTCS_2D(grid, bc, timestep);
    }

    // checking the boundary condition type first does not work
    // if the boundary condition types change dynamically for a grid
    // meaning:
    // - boundary condition type needs to be checked for every single boundary cell
    //      -> this check is last in order
    // - 
}
