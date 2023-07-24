#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <Eigen/Core>
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

/*******************************************
class WallElement {
    public:
        WallElement() {
            this->type = BC_TYPE_CLOSED;
            this->value = 0;
        }

        WallElement(double value) {
            this->type = BC_TYPE_CONSTANT;
            this->value = value;
        }

        BC_TYPE getType() {
            return this->type;
        }

        double getValue() {
            return this->value;
        }

    private:
        BC_TYPE type;
        double value;
};

class BoundaryWall {
    public: 
        BoundaryWall(int length) {
            // create array with length many wall elements
        }

    private:
        BC_SIDE side;
        int length;
        vector<WallElement> wall;
};
***********************/

class Boundary {
    public:

        /**
         * @brief Construct a new Boundary object
         * 
         * @param grid
         * @param type  
         */
        Boundary(Grid grid, BC_TYPE type);

        /**
         * @brief Get the Boundary Condition Type object
         * 
         * @return auto 
         */
        BC_TYPE getBoundaryConditionType();

        /**
         * @brief Set the Boundary Condition Value object
         * 
         * @param side 
         * @param values 
         */
        void setBoundaryConditionValue(BC_SIDE side, VectorXd values);

        /**
         * @brief Get the Boundary Condition Value object
         * 
         * @param side 
         * @return auto 
         */
        VectorXd getBoundaryConditionValue(BC_SIDE side);


    private:
        Grid grid;

        // need a way to save the bc type and value for each single 'boundary cell'
        // perhaps an array for each side with structs containing the bc type as well as a value
        // or another object that contains one boundary side 
        BC_TYPE type;
        VectorXd left, right, top, bottom;
};

#endif
