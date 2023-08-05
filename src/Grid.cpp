#include "TugUtils.hpp"
#include <tug/Grid.hpp>
#include <iostream>

Grid::Grid(int length) {
    if (length <= 3) {
        throw_invalid_argument("Given grid length too small. Must be greater than 3.");
    }

    this->row = 1;
    this->col = length;
    this->domainCol = length; // default: same size as length
    this->deltaCol = double(this->domainCol)/double(this->col); // -> 1
    this->dim = 1;

    // TODO move to the case when Simulation is set to constant and use as default
    this->concentrations = MatrixXd::Constant(1, col, 20);
    this->alphaX = MatrixXd::Constant(1, col, 1);
}

Grid::Grid(int row, int col) {
    if (row <= 3 || col <= 3) {
        throw_invalid_argument("Given grid dimensions too small. Must each be greater than 3.");
    }

    this->row = row;
    this->col = col;
    this->domainRow = row; // default: same size as row
    this->domainCol = col; // default: same size as col
    this->deltaRow = double(this->domainRow)/double(this->row); // -> 1
    this->deltaCol = double(this->domainCol)/double(this->col); // -> 1
    this->dim = 2;

    // TODO move to the case when Simulation is set to constant and use as default
    this->concentrations = MatrixXd::Constant(row, col, 20);
    this->alphaX = MatrixXd::Constant(row, col, 1);
    this->alphaY = MatrixXd::Constant(row, col, 1);
}

void Grid::setConcentrations(MatrixXd concentrations) {
    if (concentrations.rows() != this->row || concentrations.cols() != this->col) {
        throw_invalid_argument("Given matrix of concentrations mismatch with Grid dimensions!");
    }

    this->concentrations = concentrations;
}

MatrixXd Grid::getConcentrations() {
    return this->concentrations;
}

void Grid::setAlpha(MatrixXd alpha) {
    if (dim != 1) {
        throw_invalid_argument("Grid is not one dimensional, you should probably use 2D setter function!");
    }
    if (alpha.rows() != 1 || alpha.cols() != this->col) {
        throw_invalid_argument("Given matrix of alpha coefficients mismatch with Grid dimensions!");
    }

    this->alphaX = alpha;
}

void Grid::setAlpha(MatrixXd alphaX, MatrixXd alphaY) {
    if (dim != 2) {
        throw_invalid_argument("Grid is not two dimensional, you should probably use 1D setter function!");
    }
    if (alphaX.rows() != this->row || alphaX.cols() != this->col) {
        throw_invalid_argument("Given matrix of alpha coefficients in x-direction mismatch with GRid dimensions!");
    }
    if (alphaY.rows() != this->row || alphaY.cols() != this->col) {
        throw_invalid_argument("Given matrix of alpha coefficients in y-direction mismatch with GRid dimensions!");
    }

    this->alphaX = alphaX;
    this->alphaY = alphaY;
}

MatrixXd Grid::getAlpha() {
    if (dim != 1) {
        throw_invalid_argument("Grid is not one dimensional, you should probably use either getAlphaX() or getAlphaY()!");
    }

    return this->alphaX;
}

MatrixXd Grid::getAlphaX() {
    if (dim != 2) {
        throw_invalid_argument("Grid is not two dimensional, you should probably use getAlpha()!");
    }

    return this->alphaX;
}

MatrixXd Grid::getAlphaY() {
    if (dim != 2) {
        throw_invalid_argument("Grid is not two dimensional, you should probably use getAlpha()!");
    }

    return this->alphaY;
}

int Grid::getDim() {
    return dim;
}

int Grid::getLength() {
    if (dim != 1) {
        throw_invalid_argument("Grid is not one dimensional, you should probably use getRow() or getCol()!");
    }

    return col;
}

int Grid::getRow() {
    return row;
}

int Grid::getCol() {
    return col;
}

void Grid::setDomain(int domainLength) {
    if (dim != 1) {
        throw_invalid_argument("Grid is not one dimensional, you should probaly use the 2D domain setter!");
    }
    if (domainLength <= 0) {
        throw_invalid_argument("Given domain length is not positive!");
    }

    this->domainCol = domainLength;
    this->deltaCol = double(this->domainCol)/double(this->col);
}

void Grid::setDomain(int domainRow, int domainCol) {
    if (dim != 2) {
        throw_invalid_argument("Grid is not two dimensional, you should probably use the 1D domain setter!");
    }
    if (domainRow <= 0 || domainCol <= 0) {
        throw_invalid_argument("Given domain size is not positive!");
    }

    this->domainRow = domainRow;
    this->domainCol = domainCol;
    this->deltaRow = double(this->domainRow)/double(this->row);
    this->deltaCol = double(this->domainCol)/double(this->col);
}

double Grid::getDeltaCol() {
    return this->deltaCol;
}

double Grid::getDeltaRow() {
    if (dim != 2) {
        throw_invalid_argument("Grid is not two dimensional, meaning there is no delta in y-direction!");
    }

    return this->deltaRow;
}
