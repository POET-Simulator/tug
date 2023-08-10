#include "FTCS.cpp"
#include <arm_neon.h>
#include <cstddef>
#include <tug/Boundary.hpp>

using namespace Eigen;


static SparseMatrix<double> createCoeffMatrix(MatrixXd &alpha, int dim, int rowIndex, double sx) {

    // square matrix of column^2 dimension for the coefficients
    SparseMatrix<double> cm(dim, dim);
    cm.reserve(VectorXi::Constant(dim, 3));

    // left column
    cm.insert(0,0) = 1 + sx * (calcAlphaIntercell(alpha(rowIndex,0), alpha(rowIndex,1))
                                + 2 * alpha(rowIndex,0));
    cm.insert(0,1) = -sx * calcAlphaIntercell(alpha(rowIndex,0), alpha(rowIndex,1));

    // inner columns
    int n = dim-1;
    for (int i = 1; i < n; i++) {
        cm.insert(i,i-1) = -sx * calcAlphaIntercell(alpha(rowIndex,i-1), alpha(rowIndex,i));
        cm.insert(i,i) = 1 + sx * (
                            calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex,i+1))
                            + calcAlphaIntercell(alpha(rowIndex,i-1), alpha(i,1))
                            );
        cm.insert(i,i+1) = -sx * calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex,i+1));
    }

    // right column
    cm.insert(n,n-1) = -sx * calcAlphaIntercell(alpha(rowIndex,n-1), alpha(rowIndex,n));
    cm.insert(n,n) = 1 + sx * (calcAlphaIntercell(alpha(rowIndex,n-1), alpha(rowIndex,n))
                                + 2 * alpha(rowIndex,n));

    cm.makeCompressed();

    return cm;
}


static VectorXd createSolutionVector(MatrixXd &concentrations, MatrixXd &alphaX, MatrixXd &alphaY, 
                                        VectorXd &bcLeft, VectorXd &bcRight, VectorXd &bcTop, 
                                        VectorXd &bcBottom, int length, int rowIndex, double sx, double sy) {

    VectorXd sv(length);
    int numRows = concentrations.rows();

    // inner rows
    if (rowIndex > 0 && rowIndex < numRows-1) {
        for (int i = 0; i < length; i++) {
            sv(i) = sy * calcAlphaIntercell(alphaY(rowIndex,i), alphaY(rowIndex+1,i))
                            * concentrations(rowIndex+1,i)
                        + (
                            1 - sy * (
                                calcAlphaIntercell(alphaY(rowIndex,i), alphaY(rowIndex+1,i))
                                + calcAlphaIntercell(alphaY(rowIndex-1,i), alphaY(rowIndex,i))
                            )
                        ) * concentrations(rowIndex,i)
                        + sy * calcAlphaIntercell(alphaY(rowIndex-1,i), alphaY(rowIndex,i)) 
                            * concentrations(rowIndex-1,i)
                    ;
        }
    }

    // first row
    if (rowIndex == 0) {
        for (int i = 0; i < length; i++) {
            sv(i) = sy * calcAlphaIntercell(alphaY(rowIndex,i), alphaY(rowIndex+1,i))
                            * concentrations(rowIndex,i)
                        + (
                            1 - sy * (
                                calcAlphaIntercell(alphaY(rowIndex,i), alphaY(rowIndex+1,i))
                                + 2 * alphaY(rowIndex,i)
                            )
                        ) * concentrations(rowIndex,i)
                        + sy * alphaY(rowIndex,i) * bcTop(i)
                    ;
        }
    }

    // last row
    if (rowIndex == numRows-1) {
        for (int i = 0; i < length; i++) {
            sv(i) = sy * alphaY(rowIndex,i) * bcBottom(i)
                        + (
                            1 - sy * (
                                2 * alphaY(rowIndex,i)
                                + calcAlphaIntercell(alphaY(rowIndex-1,i), alphaY(rowIndex,i))
                            )
                        ) * concentrations(rowIndex,i)
                        + sy * calcAlphaIntercell(alphaY(rowIndex-1,i), alphaY(rowIndex,i)) 
                            * concentrations(rowIndex-1,i)
                    ;
        }
    }

    // first column -> additional fixed concentration change from perpendicular dimension
    sv(0) += 2 * sx * alphaX(rowIndex,0) * bcLeft(0);

    // last column -> additional fixed concentration change from perpendicular dimension
    sv(length-1) += 2 * sx * alphaX(rowIndex,length-1) * bcRight(length-1);

    return sv;
}


// BTCS solution for 1D grid
static void BTCS_1D(Grid &grid, Boundary &bc, double &timestep) {
    // TODO
}


static VectorXd solve(SparseMatrix<double> &A, VectorXd &b) {
    SparseLU<SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    return solver.solve(b);
}


// BTCS solution for 2D grid
static void BTCS_2D(Grid &grid, Boundary &bc, double &timestep) {
    int rowMax = grid.getRow();
    int colMax = grid.getCol();
    double sx = timestep / (2 * grid.getDeltaCol() * grid.getDeltaCol());
    double sy = timestep / (2 * grid.getDeltaRow() * grid.getDeltaRow());

    MatrixXd concentrations_t1 = MatrixXd::Constant(rowMax, colMax, 0);
    VectorXd row_t1(colMax);

    SparseMatrix<double> A;
    VectorXd b;

    MatrixXd alphaX = grid.getAlphaX();
    MatrixXd alphaY = grid.getAlphaY();
    VectorXd bcLeft = bc.getBoundarySideValues(BC_SIDE_LEFT);
    VectorXd bcRight = bc.getBoundarySideValues(BC_SIDE_RIGHT);
    VectorXd bcTop = bc.getBoundarySideValues(BC_SIDE_TOP);
    VectorXd bcBottom = bc.getBoundarySideValues(BC_SIDE_BOTTOM);

    MatrixXd concentrations = grid.getConcentrations();
    for (int i = 0; i < rowMax; i++) {
      
        A = createCoeffMatrix(alphaX, colMax, i, sx);
        b = createSolutionVector(concentrations, alphaX, alphaY, bcLeft, bcRight, 
                                    bcTop, bcBottom, colMax, i, sx, sy);
        row_t1 = solve(A, b);
        
        for (int j = 0; j < colMax; j++) {
            concentrations_t1(i,j) = row_t1(j);
        }
        
    }
    concentrations_t1.transposeInPlace();
    concentrations.transposeInPlace();
    alphaX.transposeInPlace();
    alphaY.transposeInPlace();
    for (int i = 0; i < colMax; i++) {

        A = createCoeffMatrix(alphaY, rowMax, i, sy);
        b = createSolutionVector(concentrations_t1, alphaY, alphaX, bcTop, bcBottom, 
                                    bcLeft, bcRight, rowMax, i, sy, sx);
        row_t1 = solve(A, b);

        for (int j = 0; j < rowMax; j++) {
            concentrations(i,j) = row_t1(j);
        }
    }
    concentrations.transposeInPlace();

    grid.setConcentrations(concentrations);
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