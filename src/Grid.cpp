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

auto Grid::getConcentrations() {
    return this->concentrations;
}

void Grid::setAlpha(Matrix2d alpha) {
    this->alpha_x = alpha;
}

void Grid::setAlpha(Matrix2d alpha_x, Matrix2d alpha_y) {
    this->alpha_x = alpha_x;
    this->alpha_y = alpha_y;
}

auto Grid::getDim() {
    return dim;
}

auto Grid::getRow() {
    return row;
}

auto Grid::getCol() {
    return col;
}