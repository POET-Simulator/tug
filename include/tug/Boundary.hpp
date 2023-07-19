#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <Eigen/Core>
#include "Grid.hpp"

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
        void setBoundaryConditionValue(BC_SIDE side, VectorXd &values);

        /**
         * @brief Get the Boundary Condition Value object
         * 
         * @param side 
         * @return auto 
         */
        VectorXd getBoundaryConditionValue(BC_SIDE side);


    private:
        Grid grid;
        BC_TYPE type;
        VectorXd left, right, top, bottom;
};

#endif
