#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <Eigen/Core>
#include <cstddef>
#include "Grid.hpp"

using namespace std;
using namespace Eigen;

enum BC_TYPE {
    BC_TYPE_CLOSED,
    BC_TYPE_CONSTANT
};

enum BC_SIDE {
    BC_SIDE_LEFT,
    BC_SIDE_RIGHT,
    BC_SIDE_TOP,
    BC_SIDE_BOTTOM
};

class BoundaryElement {
    public:
        // bc type closed
        BoundaryElement();

        // bc type constant
        BoundaryElement(double value);

        void setType(BC_TYPE type);

        void setValue(double value);

        BC_TYPE getType();

        double getValue();

    private:
        BC_TYPE type;
        double value;
};

class BoundaryWall {
    public: 
        BoundaryWall(int length);

        void setWall(BC_TYPE type, double value = NAN);

        vector<BoundaryElement> getWall();

        void setBoundaryElement(int index, BC_TYPE type, double value = NAN);

        BoundaryElement getBoundaryElement();

    private:
        BC_SIDE side;
        int length;
        vector<BoundaryElement> wall;

};

class Boundary {
    public:

        /**
         * @brief Construct a new Boundary object
         * 
         * @param grid 
         */
        Boundary(Grid grid);

        void setBoundarySideClosed(BC_SIDE side);

        void setBoundarySideConstant(BC_SIDE side, double value);

        void setBoundaryElementClosed(BC_SIDE side, int index);

        void setBoundaryElementConstant(BC_SIDE side, int index, double value);

        vector<BoundaryElement> getBoundarySide(BC_SIDE side);

        BoundaryElement getBoundaryElement(BC_SIDE side, int index);

        BC_TYPE getBoundaryElementType(BC_SIDE side, int index);

        double getBoundaryElementValue(BC_SIDE side, int index);

    private:
        Grid grid;
        
        vector<vector<BoundaryElement>> boundaries;
};

#endif

