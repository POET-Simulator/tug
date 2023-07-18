#include <tug/Grid.hpp>

Grid::Grid(int col) {
    this->col = col;

    this->dim = 1;
    this->concentrations = Matrix2d::Constant(1, col, 1);
    this->alpha_x = Matrix2d::Constant(1, col, 1);
}

Grid::Grid(int row, int col) {
    this->row = row;
    this->col = col;

    this->dim = 2;
    this->concentrations = Matrix2d::Constant(row, col, 1);
    this->alpha_x = Matrix2d::Constant(row, col, 1);
    this->alpha_y = Matrix2d::Constant(row, col, 1);
}

void Grid::setConcentrations(Matrix2d concentrations) {
    this->concentrations = concentrations;
}

Matrix2d Grid::getConcentrations() {
    return this->concentrations;
}

void Grid::setAlpha(Matrix2d alpha) {
    this->alpha_x = alpha;
}

void Grid::setAlpha(Matrix2d alpha_x, Matrix2d alpha_y) {
    this->alpha_x = alpha_x;
    this->alpha_y = alpha_y;
}

Matrix2d Grid::getAlphaX() {
    return this->alpha_x;
}

Matrix2d Grid::getAlphaY() {
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
    this->delta_col = this->domain_col/this->col;
}

void Grid::setDomain(int domain_row, int domain_col) {
    this->domain_row = domain_row;
    this->domain_col = domain_col;

    this->domain_row = this->domain_row/this->row;
    this->domain_col = this->domain_col/this->col;
}

double Grid::getDeltaCol() {
    return this->delta_col;
}

double Grid::getDeltaRow() {
    return this->delta_row;
}
