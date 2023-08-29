/**
 * @file BTCSv2.cpp
 * @brief Implementation of heterogenous BTCS (backward time-centered space) solution
 *        of diffusion equation in 1D and 2D space. Internally the alternating-direction
 *        implicit (ADI) method is used. Version 2, because Version 1 was an 
 *        implementation for the homogeneous BTCS solution. 
 * 
 */

#include "FTCS.cpp"
#include <tug/Boundary.hpp>
#include <omp.h>

#define NUM_THREADS_BTCS 10

using namespace Eigen;


// calculates coefficient for left boundary in constant case
static tuple<double, double> calcLeftBoundaryCoeffConstant(MatrixXd &alpha, int rowIndex, double sx) {
    double centerCoeff;
    double rightCoeff;

    centerCoeff = 1 + sx * (calcAlphaIntercell(alpha(rowIndex,0), alpha(rowIndex,1))
                                + 2 * alpha(rowIndex,0));
    rightCoeff = -sx * calcAlphaIntercell(alpha(rowIndex,0), alpha(rowIndex,1));

    return {centerCoeff, rightCoeff};
}


// calculates coefficient for left boundary in closed case
static tuple<double, double> calcLeftBoundaryCoeffClosed(MatrixXd &alpha, int rowIndex, double sx) {
    double centerCoeff;
    double rightCoeff;

    centerCoeff = 1 + sx * calcAlphaIntercell(alpha(rowIndex,0), alpha(rowIndex,1));
    rightCoeff = -sx * calcAlphaIntercell(alpha(rowIndex,0), alpha(rowIndex,1));

    return {centerCoeff, rightCoeff};
}


// calculates coefficient for right boundary in constant case
static tuple<double, double> calcRightBoundaryCoeffConstant(MatrixXd &alpha, int rowIndex, int n, double sx) {
    double leftCoeff;
    double centerCoeff;

    leftCoeff = -sx * calcAlphaIntercell(alpha(rowIndex,n-1), alpha(rowIndex,n));
    centerCoeff = 1 + sx * (calcAlphaIntercell(alpha(rowIndex,n-1), alpha(rowIndex,n))
                                + 2 * alpha(rowIndex,n));

    return {leftCoeff, centerCoeff};
}


// calculates coefficient for right boundary in closed case
static tuple<double, double> calcRightBoundaryCoeffClosed(MatrixXd &alpha, int rowIndex, int n, double sx) {
    double leftCoeff;
    double centerCoeff;

    leftCoeff = -sx * calcAlphaIntercell(alpha(rowIndex,n-1), alpha(rowIndex,n));
    centerCoeff = 1 + sx * calcAlphaIntercell(alpha(rowIndex,n-1), alpha(rowIndex,n));

    return {leftCoeff, centerCoeff};
}


// creates coefficient matrix for next time step from alphas in x-direction
static SparseMatrix<double> createCoeffMatrix(MatrixXd &alpha, vector<BoundaryElement> &bcLeft, vector<BoundaryElement> &bcRight, int numCols, int rowIndex, double sx) {

    // square matrix of column^2 dimension for the coefficients
    SparseMatrix<double> cm(numCols, numCols);
    cm.reserve(VectorXi::Constant(numCols, 3));

    // left column
    BC_TYPE type = bcLeft[rowIndex].getType();
    if (type == BC_TYPE_CONSTANT) {
        auto [centerCoeffTop, rightCoeffTop] = calcLeftBoundaryCoeffConstant(alpha, rowIndex, sx);
        cm.insert(0,0) = centerCoeffTop;
        cm.insert(0,1) = rightCoeffTop;
    } else if (type == BC_TYPE_CLOSED) { 
        auto [centerCoeffTop, rightCoeffTop] = calcLeftBoundaryCoeffClosed(alpha, rowIndex, sx);
        cm.insert(0,0) = centerCoeffTop;
        cm.insert(0,1) = rightCoeffTop;
    } else {
        throw_invalid_argument("Undefined Boundary Condition Type somewhere on Left or Top!");
    }

    // inner columns
    int n = numCols-1;
    for (int i = 1; i < n; i++) {
        cm.insert(i,i-1) = -sx * calcAlphaIntercell(alpha(rowIndex,i-1), alpha(rowIndex,i));
        cm.insert(i,i) = 1 + sx * (
                            calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex,i+1))
                            + calcAlphaIntercell(alpha(rowIndex,i-1), alpha(rowIndex,i))
                            )
                        ;
        cm.insert(i,i+1) = -sx * calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex,i+1));
    }

    // right column
    type = bcRight[rowIndex].getType();
    if (type == BC_TYPE_CONSTANT) {
        auto [leftCoeffBottom, centerCoeffBottom] = calcRightBoundaryCoeffConstant(alpha, rowIndex, n, sx);
        cm.insert(n,n-1) = leftCoeffBottom;
        cm.insert(n,n) = centerCoeffBottom;
    } else if (type == BC_TYPE_CLOSED) { 
        auto [leftCoeffBottom, centerCoeffBottom] = calcRightBoundaryCoeffClosed(alpha, rowIndex, n, sx);
        cm.insert(n,n-1) = leftCoeffBottom;
        cm.insert(n,n) = centerCoeffBottom;
    } else {
        throw_invalid_argument("Undefined Boundary Condition Type somewhere on Right or Bottom!");
    }

    cm.makeCompressed(); // important for Eigen solver

    return cm;
}


// calculates explicity concentration at top boundary in constant case
static double calcExplicitConcentrationsTopBoundaryConstant(MatrixXd &concentrations, 
                            MatrixXd &alpha, vector<BoundaryElement> &bcTop, int rowIndex, int i, double sy) {
    double c;

    c = sy * calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex+1,i))
                            * concentrations(rowIndex,i)
                        + (
                            1 - sy * (
                                calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex+1,i))
                                + 2 * alpha(rowIndex,i)
                            )
                        ) * concentrations(rowIndex,i)
                        + sy * alpha(rowIndex,i) * bcTop[i].getValue();

    return c;
}


// calculates explicit concentration at top boundary in closed case
static double calcExplicitConcentrationsTopBoundaryClosed(MatrixXd &concentrations, 
                            MatrixXd &alpha, int rowIndex, int i, double sy) {
    double c;

    c = sy * calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex+1,i))
                            * concentrations(rowIndex,i)
                        + (
                            1 - sy * (
                                calcAlphaIntercell(alpha(rowIndex,i), alpha(rowIndex+1,i))
                            )
                        ) * concentrations(rowIndex,i);

    return c;
}


// calculates explicit concentration at bottom boundary in constant case
static double calcExplicitConcentrationsBottomBoundaryConstant(MatrixXd &concentrations, 
                            MatrixXd &alpha, vector<BoundaryElement> &bcBottom, int rowIndex, int i, double sy) {
    double c;

    c = sy * alpha(rowIndex,i) * bcBottom[i].getValue()
                        + (
                            1 - sy * (
                                2 * alpha(rowIndex,i)
                                + calcAlphaIntercell(alpha(rowIndex-1,i), alpha(rowIndex,i))
                            )
                        ) * concentrations(rowIndex,i)
                        + sy * calcAlphaIntercell(alpha(rowIndex-1,i), alpha(rowIndex,i)) 
                            * concentrations(rowIndex-1,i);

    return c;
}


// calculates explicit concentration at bottom boundary in closed case
static double calcExplicitConcentrationsBottomBoundaryClosed(MatrixXd &concentrations, 
                            MatrixXd &alpha, int rowIndex, int i, double sy) {
    double c;

    c = (
            1 - sy * (
                + calcAlphaIntercell(alpha(rowIndex-1,i), alpha(rowIndex,i))
            )
        ) * concentrations(rowIndex,i)
        + sy * calcAlphaIntercell(alpha(rowIndex-1,i), alpha(rowIndex,i)) 
            * concentrations(rowIndex-1,i);

    return c;
}


// creates a solution vector for next time step from the current state of concentrations
static VectorXd createSolutionVector(MatrixXd &concentrations, MatrixXd &alphaX, MatrixXd &alphaY, 
                                        vector<BoundaryElement> &bcLeft, vector<BoundaryElement> &bcRight, 
                                        vector<BoundaryElement> &bcTop, vector<BoundaryElement> &bcBottom, 
                                        int length, int rowIndex, double sx, double sy) {

    VectorXd sv(length);
    int numRows = concentrations.rows();
    BC_TYPE type;

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
    else if (rowIndex == 0) {
        for (int i = 0; i < length; i++) {
            type = bcTop[i].getType();
            if (type == BC_TYPE_CONSTANT) {
                sv(i) = calcExplicitConcentrationsTopBoundaryConstant(concentrations, alphaY, bcTop, rowIndex, i, sy);
            } else if (type == BC_TYPE_CLOSED) {
                sv(i) = calcExplicitConcentrationsTopBoundaryClosed(concentrations, alphaY, rowIndex, i, sy);
            } else {
                throw_invalid_argument("Undefined Boundary Condition Type somewhere on Left or Top!");
            }
        }
    }

    // last row
    else if (rowIndex == numRows-1) {
        for (int i = 0; i < length; i++) {
            type = bcBottom[i].getType();
            if (type == BC_TYPE_CONSTANT) {
                sv(i) = calcExplicitConcentrationsBottomBoundaryConstant(concentrations, alphaY, bcBottom, rowIndex, i, sy);
            } else if (type == BC_TYPE_CLOSED) {
                sv(i) = calcExplicitConcentrationsBottomBoundaryClosed(concentrations, alphaY, rowIndex, i, sy);
            } else {
                throw_invalid_argument("Undefined Boundary Condition Type somewhere on Right or Bottom!");
            }
        }
    }

    // first column -> additional fixed concentration change from perpendicular dimension in constant bc case
    if (bcLeft[rowIndex].getType() == BC_TYPE_CONSTANT) {
        sv(0) += 2 * sx * alphaX(rowIndex,0) * bcLeft[rowIndex].getValue();
    }

    // last column -> additional fixed concentration change from perpendicular dimension in constant bc case
    if (bcRight[rowIndex].getType() == BC_TYPE_CONSTANT) {
        sv(length-1) += 2 * sx * alphaX(rowIndex,length-1) * bcRight[rowIndex].getValue();
    }

    return sv;
}


// solver for linear equation system; A corresponds to coefficient matrix,
// b to the solution vector
// use of EigenLU solver
static VectorXd EigenLUAlgorithm(SparseMatrix<double> &A, VectorXd &b) {

    SparseLU<SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    return solver.solve(b);
}

// solver for linear equation system; A corresponds to coefficient matrix, 
// b to the solution vector
// implementation of Thomas Algorithm
static VectorXd ThomasAlgorithm(SparseMatrix<double> &A, VectorXd &b) {
    uint32_t n = b.size();

    Eigen::VectorXd a_diag(n);
    Eigen::VectorXd b_diag(n);
    Eigen::VectorXd c_diag(n);
    Eigen::VectorXd x_vec = b;

    // Fill diagonals vectors
    b_diag[0] = A.coeff(0, 0);
    c_diag[0] = A.coeff(0, 1);

    for (int i = 1; i < n - 1; i++) {
        a_diag[i] = A.coeff(i, i - 1);
        b_diag[i] = A.coeff(i, i);
        c_diag[i] = A.coeff(i, i + 1);
    }
    a_diag[n - 1] = A.coeff(n - 1, n - 2);
    b_diag[n - 1] = A.coeff(n - 1, n - 1);

    // start solving - c_diag and x_vec are overwritten
    n--;
    c_diag[0] /= b_diag[0];
    x_vec[0] /= b_diag[0];

    for (int i = 1; i < n; i++) {
        c_diag[i] /= b_diag[i] - a_diag[i] * c_diag[i - 1];
        x_vec[i] = (x_vec[i] - a_diag[i] * x_vec[i - 1]) /
                    (b_diag[i] - a_diag[i] * c_diag[i - 1]);
    }

    x_vec[n] = (x_vec[n] - a_diag[n] * x_vec[n - 1]) /
                (b_diag[n] - a_diag[n] * c_diag[n - 1]);

    for (int i = n; i-- > 0;) {
        x_vec[i] -= c_diag[i] * x_vec[i + 1];
    }

    return x_vec;
}


// BTCS solution for 1D grid 
static void BTCS_1D(Grid &grid, Boundary &bc, double timestep, VectorXd (*solverFunc) (SparseMatrix<double> &A, VectorXd &b)) {
    int length = grid.getLength();
    double sx = timestep / (grid.getDelta() * grid.getDelta());

    VectorXd concentrations_t1(length);

    SparseMatrix<double> A;
    VectorXd b(length);

    MatrixXd alpha = grid.getAlpha();
    vector<BoundaryElement> bcLeft = bc.getBoundarySide(BC_SIDE_LEFT);
    vector<BoundaryElement> bcRight = bc.getBoundarySide(BC_SIDE_RIGHT);

    MatrixXd concentrations = grid.getConcentrations();
    int rowIndex = 0;
    A = createCoeffMatrix(alpha, bcLeft, bcRight, length, rowIndex, sx); // this is exactly same as in 2D
    for (int i = 0; i < length; i++) {
        b(i) = concentrations(0,i);
    }
    if (bc.getBoundaryElementType(BC_SIDE_LEFT, 0) == BC_TYPE_CONSTANT) {
        b(0) += 2 * sx * alpha(0,0) * bcLeft[0].getValue();
    }
    if (bc.getBoundaryElementType(BC_SIDE_RIGHT, 0) == BC_TYPE_CONSTANT) {
        b(length-1) += 2 * sx * alpha(0,length-1) * bcRight[0].getValue();
    }

    concentrations_t1 = solverFunc(A, b);

    for (int j = 0; j < length; j++) {
        concentrations(0,j) = concentrations_t1(j);
    }

    grid.setConcentrations(concentrations);
}


// BTCS solution for 2D grid
static void BTCS_2D(Grid &grid, Boundary &bc, double timestep, VectorXd (*solverFunc) (SparseMatrix<double> &A, VectorXd &b), int numThreads) {
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
    vector<BoundaryElement> bcLeft = bc.getBoundarySide(BC_SIDE_LEFT);
    vector<BoundaryElement> bcRight = bc.getBoundarySide(BC_SIDE_RIGHT);
    vector<BoundaryElement> bcTop = bc.getBoundarySide(BC_SIDE_TOP);
    vector<BoundaryElement> bcBottom = bc.getBoundarySide(BC_SIDE_BOTTOM);

    MatrixXd concentrations = grid.getConcentrations();

    #pragma omp parallel for num_threads(numThreads) private(A, b, row_t1)
    for (int i = 0; i < rowMax; i++) {
      
        
        A = createCoeffMatrix(alphaX, bcLeft, bcRight, colMax, i, sx);
        b = createSolutionVector(concentrations, alphaX, alphaY, bcLeft, bcRight, 
                                    bcTop, bcBottom, colMax, i, sx, sy);
        
        SparseLU<SparseMatrix<double>> solver; 

        row_t1 = solverFunc(A, b);
        
        concentrations_t1.row(i) = row_t1;
    }
    
    concentrations_t1.transposeInPlace();
    concentrations.transposeInPlace();
    alphaX.transposeInPlace();
    alphaY.transposeInPlace();
    
    #pragma omp parallel for num_threads(numThreads) private(A, b, row_t1)
    for (int i = 0; i < colMax; i++) {
        // swap alphas, boundary conditions and sx/sy for column-wise calculation
        A = createCoeffMatrix(alphaY, bcTop, bcBottom, rowMax, i, sy);
        b = createSolutionVector(concentrations_t1, alphaY, alphaX, bcTop, bcBottom, 
                                    bcLeft, bcRight, rowMax, i, sy, sx);
        
        row_t1 = solverFunc(A, b);

        concentrations.row(i) = row_t1;   
    }
    
    concentrations.transposeInPlace();

    grid.setConcentrations(concentrations);
}


// entry point for EigenLU solver; differentiate between 1D and 2D grid
static void BTCS_LU(Grid &grid, Boundary &bc, double timestep, int numThreads) {
    if (grid.getDim() == 1) {
        BTCS_1D(grid, bc, timestep, EigenLUAlgorithm);
    } else if (grid.getDim() == 2) {
        BTCS_2D(grid, bc, timestep, EigenLUAlgorithm, numThreads);
    } else {
        throw_invalid_argument("Error: Only 1- and 2-dimensional grids are defined!");
    }
}

// entry point for Thomas algorithm solver; differentiate 1D and 2D grid
static void BTCS_Thomas(Grid &grid, Boundary &bc, double timestep, int numThreads) {
    if (grid.getDim() == 1) {
        BTCS_1D(grid, bc, timestep, ThomasAlgorithm);
    } else if (grid.getDim() == 2) {
        BTCS_2D(grid, bc, timestep, ThomasAlgorithm, numThreads);
    } else {
        throw_invalid_argument("Error: Only 1- and 2-dimensional grids are defined!");
    }
}