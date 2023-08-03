#include "TugUtils.hpp"
#include "tug/BoundaryCondition.hpp"
#include <iostream>
#include <omp.h>
#include <tug/Boundary.hpp>
#include <stdexcept>

using namespace std;

BoundaryElement::BoundaryElement() {
    
    this->type = BC_TYPE_CLOSED;
    this->value = NAN;
}

BoundaryElement::BoundaryElement(double value) {
    this->type = BC_TYPE_CONSTANT;
    this->value = value;
}

void BoundaryElement::setType(BC_TYPE type) {
    this->type = type;
}

void BoundaryElement::setValue(double value) {
    this->value = value;
}

BC_TYPE BoundaryElement::getType() {
    return this->type;
}

double BoundaryElement::getValue() {
    return this->value;
}

Boundary::Boundary(Grid grid) : grid(grid) {
    //probably to DEBUG assignment grid

    if (grid.getDim() == 1) {
        this->boundaries = vector<vector<BoundaryElement>>(2);

        this->boundaries[BC_SIDE_LEFT].push_back(BoundaryElement());
        this->boundaries[BC_SIDE_RIGHT].push_back(BoundaryElement());
    } else if (grid.getDim() == 2) {
        this->boundaries = vector<vector<BoundaryElement>>(4);

        this->boundaries[BC_SIDE_LEFT] = vector<BoundaryElement>(grid.getRow(), BoundaryElement());
        this->boundaries[BC_SIDE_RIGHT] = vector<BoundaryElement>(grid.getRow(), BoundaryElement());
        this->boundaries[BC_SIDE_TOP] = vector<BoundaryElement>(grid.getCol(), BoundaryElement());
        this->boundaries[BC_SIDE_BOTTOM] = vector<BoundaryElement>(grid.getCol(), BoundaryElement());
    }
}

void Boundary::setBoundarySideClosed(BC_SIDE side) {
    this->boundaries[side] = vector<BoundaryElement>(grid.getRow(), BoundaryElement());
}

void Boundary::setBoundarySideConstant(BC_SIDE side, double value) {
    this->boundaries[side] = vector<BoundaryElement>(grid.getRow(), BoundaryElement(value));
}

void Boundary::setBoundaryElementClosed(BC_SIDE side, int index) {
    // tests whether the index really points to an element of the boundary side.
    if((boundaries[side].size() < index) || index < 0){
        throw_invalid_argument("Index is selected either too large or too small.");
    }
    this->boundaries[side][index].setType(BC_TYPE_CLOSED);
}

void Boundary::setBoundaryElementConstant(BC_SIDE side, int index, double value) {
    // tests whether the index really points to an element of the boundary side.
    if((boundaries[side].size() < index) || index < 0){
        throw_invalid_argument("Index is selected either too large or too small.");
    }
    this->boundaries[side][index].setType(BC_TYPE_CONSTANT);
    this->boundaries[side][index].setValue(value);
}

vector<BoundaryElement> Boundary::getBoundarySide(BC_SIDE side) {
    return this->boundaries[side];
}

BoundaryElement Boundary::getBoundaryElement(BC_SIDE side, int index) {
    if((boundaries[side].size() < index) || index < 0){
        throw_invalid_argument("Index is selected either too large or too small.");
    }
    return this->boundaries[side][index];
}

BC_TYPE Boundary::getBoundaryElementType(BC_SIDE side, int index) {
    if((boundaries[side].size() < index) || index < 0){
        throw_invalid_argument("Index is selected either too large or too small.");
    }
    return this->boundaries[side][index].getType();
}

double Boundary::getBoundaryElementValue(BC_SIDE side, int index) {
    if((boundaries[side].size() < index) || index < 0){
        throw_invalid_argument("Index is selected either too large or too small.");
    }
    if(boundaries[side][index].getType() != BC_TYPE_CONSTANT){
        throw_invalid_argument("A value can only be output if it is a constant boundary condition.");
    }
    return this->boundaries[side][index].getValue();
}

