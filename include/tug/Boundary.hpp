#include <Eigen/Core>

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
         * @param dim 
         */
        Boundary(int dim);

        /**
         * @brief Construct a new Boundary object
         * 
         * @param dim 
         * @param type 
         */
        Boundary(int dim, BC_TYPE type);

        /**
         * @brief Set the Boundary Condition Type object
         * 
         * @param type 
         */
        void setBoundaryConditionType(BC_TYPE type);

        /**
         * @brief Get the Boundary Condition Type object
         * 
         * @return auto 
         */
        auto getBoundaryConditionType();

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
        auto getBoundaryConditionValue(BC_SIDE side);


    private:

        int dim;
        BC_TYPE type;
        VectorXd left, right, top, bottom;
};
