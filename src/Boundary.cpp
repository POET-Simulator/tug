#include "TugUtils.hpp"

#include <iostream>
#include <omp.h>
#include <stdexcept>
#include <tug/Boundary.hpp>

BoundaryElement::BoundaryElement() {

  this->type = BC_TYPE_CLOSED;
  this->value = -1; // without meaning in closed case
}

BoundaryElement::BoundaryElement(double value) {
  this->type = BC_TYPE_CONSTANT;
  this->value = value;
}

void BoundaryElement::setType(BC_TYPE type) { this->type = type; }

void BoundaryElement::setValue(double value) {
  if (value < 0) {
    throw_invalid_argument("No negative concentration allowed.");
  }
  if (type == BC_TYPE_CLOSED) {
    throw_invalid_argument(
        "No constant boundary concentrations can be set for closed "
        "boundaries. Please change type first.");
  }
  this->value = value;
}

BC_TYPE BoundaryElement::getType() { return this->type; }

double BoundaryElement::getValue() { return this->value; }

Boundary::Boundary(Grid grid) : grid(grid) {
  if (grid.getDim() == 1) {
    this->boundaries = std::vector<std::vector<BoundaryElement>>(
        2); // in 1D only left and right boundary

    this->boundaries[BC_SIDE_LEFT].push_back(BoundaryElement());
    this->boundaries[BC_SIDE_RIGHT].push_back(BoundaryElement());
  } else if (grid.getDim() == 2) {
    this->boundaries = std::vector<std::vector<BoundaryElement>>(4);

    this->boundaries[BC_SIDE_LEFT] =
        std::vector<BoundaryElement>(grid.getRow(), BoundaryElement());
    this->boundaries[BC_SIDE_RIGHT] =
        std::vector<BoundaryElement>(grid.getRow(), BoundaryElement());
    this->boundaries[BC_SIDE_TOP] =
        std::vector<BoundaryElement>(grid.getCol(), BoundaryElement());
    this->boundaries[BC_SIDE_BOTTOM] =
        std::vector<BoundaryElement>(grid.getCol(), BoundaryElement());
  }
}

void Boundary::setBoundarySideClosed(BC_SIDE side) {
  if (grid.getDim() == 1) {
    if ((side == BC_SIDE_BOTTOM) || (side == BC_SIDE_TOP)) {
      throw_invalid_argument(
          "For the one-dimensional case, only the BC_SIDE_LEFT and "
          "BC_SIDE_RIGHT borders exist.");
    }
  }

  int n;
  if (side == BC_SIDE_LEFT || side == BC_SIDE_RIGHT) {
    n = grid.getRow();
  } else {
    n = grid.getCol();
  }
  this->boundaries[side] = std::vector<BoundaryElement>(n, BoundaryElement());
}

void Boundary::setBoundarySideConstant(BC_SIDE side, double value) {
  if (grid.getDim() == 1) {
    if ((side == BC_SIDE_BOTTOM) || (side == BC_SIDE_TOP)) {
      throw_invalid_argument(
          "For the one-dimensional case, only the BC_SIDE_LEFT and "
          "BC_SIDE_RIGHT borders exist.");
    }
  }

  int n;
  if (side == BC_SIDE_LEFT || side == BC_SIDE_RIGHT) {
    n = grid.getRow();
  } else {
    n = grid.getCol();
  }
  this->boundaries[side] = std::vector<BoundaryElement>(n, BoundaryElement(value));
}

void Boundary::setBoundaryElementClosed(BC_SIDE side, int index) {
  // tests whether the index really points to an element of the boundary side.
  if ((boundaries[side].size() < index) || index < 0) {
    throw_invalid_argument("Index is selected either too large or too small.");
  }
  this->boundaries[side][index].setType(BC_TYPE_CLOSED);
}

void Boundary::setBoundaryElementConstant(BC_SIDE side, int index,
                                          double value) {
  // tests whether the index really points to an element of the boundary side.
  if ((boundaries[side].size() < index) || index < 0) {
    throw_invalid_argument("Index is selected either too large or too small.");
  }
  this->boundaries[side][index].setType(BC_TYPE_CONSTANT);
  this->boundaries[side][index].setValue(value);
}

const std::vector<BoundaryElement> Boundary::getBoundarySide(BC_SIDE side) {
  if (grid.getDim() == 1) {
    if ((side == BC_SIDE_BOTTOM) || (side == BC_SIDE_TOP)) {
      throw_invalid_argument(
          "For the one-dimensional trap, only the BC_SIDE_LEFT and "
          "BC_SIDE_RIGHT borders exist.");
    }
  }
  return this->boundaries[side];
}

Eigen::VectorXd Boundary::getBoundarySideValues(BC_SIDE side) {
  int length = boundaries[side].size();
  Eigen::VectorXd values(length);

  for (int i = 0; i < length; i++) {
    if (getBoundaryElementType(side, i) == BC_TYPE_CLOSED) {
      values(i) = -1;
      continue;
    }
    values(i) = getBoundaryElementValue(side, i);
  }

  return values;
}

BoundaryElement Boundary::getBoundaryElement(BC_SIDE side, int index) {
  if ((boundaries[side].size() < index) || index < 0) {
    throw_invalid_argument("Index is selected either too large or too small.");
  }
  return this->boundaries[side][index];
}

BC_TYPE Boundary::getBoundaryElementType(BC_SIDE side, int index) {
  if ((boundaries[side].size() < index) || index < 0) {
    throw_invalid_argument("Index is selected either too large or too small.");
  }
  return this->boundaries[side][index].getType();
}

double Boundary::getBoundaryElementValue(BC_SIDE side, int index) {
  if ((boundaries[side].size() < index) || index < 0) {
    throw_invalid_argument("Index is selected either too large or too small.");
  }
  if (boundaries[side][index].getType() != BC_TYPE_CONSTANT) {
    throw_invalid_argument(
        "A value can only be output if it is a constant boundary condition.");
  }
  return this->boundaries[side][index].getValue();
}
