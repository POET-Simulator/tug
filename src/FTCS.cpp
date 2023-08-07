/**
 * @file FTCS.cpp
 * @brief Implementation of heterogenous FTCS (forward time-centered space) solution
 *        of diffusion equation in 1D and 2D space.
 * 
 */

#include "TugUtils.hpp"
#include <cstddef>
#include <tug/Boundary.hpp>
#include <iostream>

using namespace std;


// calculates arithmetic or harmonic mean of alpha between two cells
static double calcAlphaIntercell(double &alpha1, double &alpha2, bool useHarmonic = true) {
    if (useHarmonic) {
        return double(2) / ((double(1)/alpha1) + (double(1)/alpha2));
    } else {
        return 0.5 * (alpha1 + alpha2);
    }
}


// calculates horizontal change on one cell independent of boundary type
static double calcHorizontalChange(Grid &grid, int &row, int &col) {

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


// calculates vertical change on one cell independent of boundary type
static double calcVerticalChange(Grid &grid, int &row, int &col) {
    
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


// calculates horizontal change on one cell with a constant left boundary
static double calcHorizontalChangeLeftBoundaryConstant(Grid &grid, Boundary &bc, int &row, int &col) {

    double result = 
        calcAlphaIntercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
            * grid.getConcentrations()(row,col+1)
        - (
            calcAlphaIntercell(grid.getAlphaX()(row,col+1), grid.getAlphaX()(row,col)) 
            + 2 * grid.getAlphaX()(row,col)
            ) 
            * grid.getConcentrations()(row,col) 
        + 2 * grid.getAlphaX()(row,col) * bc.getBoundaryElementValue(BC_SIDE_LEFT, row);

    return result;
}


// calculates horizontal change on one cell with a closed left boundary
static double calcHorizontalChangeLeftBoundaryClosed(Grid &grid, int &row, int &col) {
    
    double result = 
        calcAlphaIntercell(grid.getAlphaX()(row, col+1), grid.getAlphaX()(row, col)) 
            * (grid.getConcentrations()(row, col+1) - grid.getConcentrations()(row, col));
    
    return result;
}


// checks boundary condition type for a cell on the left edge of grid
static double calcHorizontalChangeLeftBoundary(Grid &grid, Boundary &bc, int &row, int &col) {
    if (bc.getBoundaryElementType(BC_SIDE_LEFT, col) == BC_TYPE_CONSTANT) {
        return calcHorizontalChangeLeftBoundaryConstant(grid, bc, row, col);
    } else if (bc.getBoundaryElementType(BC_SIDE_LEFT, col) == BC_TYPE_CLOSED) {
        return calcHorizontalChangeLeftBoundaryClosed(grid, row, col);
    } else {
        throw_invalid_argument("Undefined Boundary Condition Type!");
    }
}


// calculates horizontal change on one cell with a constant right boundary
static double calcHorizontalChangeRightBoundaryConstant(Grid &grid, Boundary &bc, int &row, int &col) {

    double result = 
        2 * grid.getAlphaX()(row,col) * bc.getBoundaryElementValue(BC_SIDE_RIGHT, row)
        - (
            calcAlphaIntercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
            + 2 * grid.getAlphaX()(row,col)
            ) 
            * grid.getConcentrations()(row,col)
        + calcAlphaIntercell(grid.getAlphaX()(row,col-1), grid.getAlphaX()(row,col))
            * grid.getConcentrations()(row,col-1);

    return result;
}


// calculates horizontal change on one cell with a closed right boundary
static double calcHorizontalChangeRightBoundaryClosed(Grid &grid, int &row, int &col) {
    
    double result = 
        - (calcAlphaIntercell(grid.getAlphaX()(row, col-1), grid.getAlphaX()(row, col))
         * (grid.getConcentrations()(row, col) - grid.getConcentrations()(row, col-1)));
    
    return result;
}


// checks boundary condition type for a cell on the right edge of grid
static double calcHorizontalChangeRightBoundary(Grid &grid, Boundary &bc, int &row, int &col) {
    if (bc.getBoundaryElementType(BC_SIDE_RIGHT, col) == BC_TYPE_CONSTANT) {
        return calcHorizontalChangeRightBoundaryConstant(grid, bc, row, col);
    } else if (bc.getBoundaryElementType(BC_SIDE_RIGHT, col) == BC_TYPE_CLOSED) {
        return calcHorizontalChangeRightBoundaryClosed(grid, row, col);
    } else {
        throw_invalid_argument("Undefined Boundary Condition Type!");
    }
}


// calculates vertical change on one cell with a constant top boundary
static double calcVerticalChangeTopBoundaryConstant(Grid &grid, Boundary &bc, int &row, int &col) {
    
    double result = 
        calcAlphaIntercell(grid.getAlphaY()(row+1, col), grid.getAlphaY()(row, col)) 
            * grid.getConcentrations()(row+1,col)
        - (
            calcAlphaIntercell(grid.getAlphaY()(row+1, col), grid.getAlphaY()(row, col)) 
            + 2 * grid.getAlphaY()(row, col)
            ) 
            * grid.getConcentrations()(row, col)
        + 2 * grid.getAlphaY()(row, col) * bc.getBoundaryElementValue(BC_SIDE_TOP, col);

    return result;
}


// calculates vertical change on one cell with a closed top boundary
static double calcVerticalChangeTopBoundaryClosed(Grid &grid, int &row, int &col) {
    
    double result = 
        calcAlphaIntercell(grid.getAlphaY()(row+1, col), grid.getConcentrations()(row, col)) 
            * (grid.getConcentrations()(row+1, col) - grid.getConcentrations()(row, col));

    return result;
}


// checks boundary condition type for a cell on the top edge of grid
static double calcVerticalChangeTopBoundary(Grid &grid, Boundary &bc, int &row, int &col) {
    if (bc.getBoundaryElementType(BC_SIDE_TOP, col) == BC_TYPE_CONSTANT) {
        return calcVerticalChangeTopBoundaryConstant(grid, bc, row, col);
    } else if (bc.getBoundaryElementType(BC_SIDE_TOP, col) == BC_TYPE_CLOSED) {
        return calcVerticalChangeTopBoundaryClosed(grid, row, col);
    } else {
        throw_invalid_argument("Undefined Boundary Condition Type!");
    }
}


// calculates vertical change on one cell with a constant bottom boundary
static double calcVerticalChangeBottomBoundaryConstant(Grid &grid, Boundary &bc, int &row, int &col) {

    double result = 
        2 * grid.getAlphaY()(row, col) * bc.getBoundaryElementValue(BC_SIDE_BOTTOM, col)
        - (
            calcAlphaIntercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)) 
            + 2 * grid.getAlphaY()(row, col)
            ) 
            * grid.getConcentrations()(row, col)
        + calcAlphaIntercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col)) 
            * grid.getConcentrations()(row-1,col);

    return result;
}


// calculates vertical change on one cell with a closed bottom boundary
static double calcVerticalChangeBottomBoundaryClosed(Grid &grid, int &row, int &col) {

    double result = 
        - (calcAlphaIntercell(grid.getAlphaY()(row, col), grid.getAlphaY()(row-1, col))
            * (grid.getConcentrations()(row, col) - grid.getConcentrations()(row-1, col)));

    return result;
}


// checks boundary condition type for a cell on the bottom edge of grid
static double calcVerticalChangeBottomBoundary(Grid &grid, Boundary &bc, int &row, int &col) {
    if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, col) == BC_TYPE_CONSTANT) {
        return calcVerticalChangeBottomBoundaryConstant(grid, bc, row, col);
    } else if (bc.getBoundaryElementType(BC_SIDE_BOTTOM, col) == BC_TYPE_CLOSED) {
        return calcVerticalChangeBottomBoundaryClosed(grid, row, col);
    } else {
        throw_invalid_argument("Undefined Boundary Condition Type!");
    }
}


// FTCS solution to 1D grid
static void FTCS_1D(Grid &grid, Boundary &bc, double &timestep) {
    int colMax = grid.getCol();
    double deltaCol = grid.getDeltaCol();

    // matrix for concentrations at time t+1
    MatrixXd concentrations_t1 = MatrixXd::Constant(1, colMax, 0);

    // only one row in 1D case -> row constant at index 0
    int row = 0;

    // inner cells
    // independent of boundary condition type
    for (int col = 1; col < colMax-1; col++) {
        concentrations_t1(row,col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol)
                * (
                    calcHorizontalChange(grid, row, col)
                )
            ;
    }

    // left boundary; hold column constant at index 0
    int col = 0;
    concentrations_t1(row, col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChangeLeftBoundary(grid, bc, row, col)
                )
            ;


    // right boundary; hold column constant at max index
    col = colMax-1;
    concentrations_t1(row,col) = grid.getConcentrations()(row,col)
            + timestep / (deltaCol*deltaCol) 
                * (
                    calcHorizontalChangeRightBoundary(grid, bc, row, col)
                )
            ;

    // overwrite obsolete concentrations
    grid.setConcentrations(concentrations_t1);
}


// FTCS solution to 2D grid
static void FTCS_2D(Grid &grid, Boundary &bc, double &timestep) {
    int rowMax = grid.getRow();
    int colMax = grid.getCol();
    double deltaRow = grid.getDeltaRow();
    double deltaCol = grid.getDeltaCol();

    // matrix for concentrations at time t+1
    MatrixXd concentrations_t1 = MatrixXd::Constant(rowMax, colMax, 0);

    // inner cells
    // these are independent of the boundary condition type
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
    // hold column constant at index 0
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

    // right without corners / looping over rows
    // hold column constant at max index
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


    // top without corners / looping over columns
    // hold row constant at index 0
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

    // bottom without corners / looping over columns
    // hold row constant at max index
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
    // hold row and column constant at 0
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
    // hold row constant at 0 and column constant at max index
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
    // hold row constant at max index and column constant at 0
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
    // hold row and column constant at max index
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

    // overwrite obsolete concentrations 
    grid.setConcentrations(concentrations_t1);
}


// entry point; differentiate between 1D and 2D grid
static void FTCS(Grid &grid, Boundary &bc, double &timestep) {
    if (grid.getDim() == 1) {
        FTCS_1D(grid, bc, timestep);
    } else {
        FTCS_2D(grid, bc, timestep);
    }
}
