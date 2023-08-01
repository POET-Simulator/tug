#include "TugUtils.hpp"
#include <tug/Grid.hpp>
#include <iostream>

Grid::Grid(int col) {
    this->col = col;
    this->domain_col = col;
    this->delta_col = double(this->domain_col)/double(this->col);

    this->dim = 1;
    this->concentrations = MatrixXd::Constant(1, col, 20);
    this->alpha_x = MatrixXd::Constant(1, col, 1);
}

Grid::Grid(int row, int col) {
    // TODO check for reasonable dimensions
    if (row < 1 || col < 1) {
        throw_invalid_argument("Either row or col too small!");
    }

    this->row = row;
    this->col = col;
    this->domain_row = row;
    this->domain_col = col;
    this->delta_row = double(this->domain_row)/double(this->row);
    this->delta_col = double(this->domain_col)/double(this->col);

    this->dim = 2;
    this->concentrations = MatrixXd::Constant(row, col, 20);
    this->alpha_x = MatrixXd::Constant(row, col, 1);
    this->alpha_y = MatrixXd::Constant(row, col, 1);
}

void Grid::setConcentrations(MatrixXd concentrations) {
    this->concentrations = concentrations;
}

MatrixXd Grid::getConcentrations() {
    return this->concentrations;
}

void Grid::setAlpha(MatrixXd alpha) {
    this->alpha_x = alpha;
}

void Grid::setAlpha(MatrixXd alpha_x, MatrixXd alpha_y) {
    this->alpha_x = alpha_x;
    this->alpha_y = alpha_y;
}

MatrixXd Grid::getAlphaX() {
    return this->alpha_x;
}

MatrixXd Grid::getAlphaY() {
    return this->alpha_y;
}

int Grid::getDim() {
    return dim;
}

int Grid::getRow() {
    return row;
}

int Grid::getCol() {
    return col;
}

void Grid::setDomain(int domain_col) {
    this->domain_col = domain_col;
    this->delta_col = double(this->domain_col)/this->col;
}

void Grid::setDomain(int domain_row, int domain_col) {
    this->domain_row = domain_row;
    this->domain_col = domain_col;

    this->delta_row = double(this->domain_row)/double(this->row);
    this->delta_col = double(this->domain_col)/double(this->col);
}

double Grid::getDeltaCol() {
    return this->delta_col;
}

double Grid::getDeltaRow() {
    return this->delta_row;
}
