/**
 * @file Grid.hpp
 * @brief API of Grid class, that holds a matrix with concenctrations and a 
 *        respective matrix/matrices of alpha coefficients. 
 * 
 */

#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace Eigen;

class Grid {
    public:

        /**
         * @brief Constructs a new 1D-Grid object of a given length, which holds a matrix 
         *        with concentrations and a respective matrix of alpha coefficients. 
         *        The domain length is per default the same as the length. The concentrations
         *        are all 20 by default and the alpha coefficients are 1. 
         * 
         * @param length Length of the 1D-Grid. Must be greater than 3. 
         */
        Grid(int length);

        /**
        * @brief Constructs a new 2D-Grid object of given dimensions, which holds a matrix
        *        with concentrations and the respective matrices of alpha coefficient for 
        *        each direction. The domain in x- and y-direction is per default equal to 
        *        the col length and row length, respectively. 
        *        The concentrations are all 20 by default across the entire grid and the 
        *        alpha coefficients 1 in both directions. 
        * 
        * @param row Length of the 2D-Grid in y-direction. Must be greater than 3. 
        * @param col Length of the 2D-Grid in x-direction. Must be greater than 3. 
        */
        Grid(int row, int col);

        /**
        * @brief Sets the concentrations matrix for a 1D or 2D-Grid.
        * 
        * @param concentrations An Eigen3 MatrixXd holding the concentrations. Matrix must 
        *                       have correct dimensions as defined in row and col. (Or length, 
        *                       in 1D case). 
        */
        void setConcentrations(MatrixXd concentrations);

        /**
         * @brief Gets the concentrations matrix for a Grid. 
         * 
         * @return MatrixXd An Eigen3 matrix holding the concentrations and having the 
         *                  same dimensions as the grid. 
         */
        const MatrixXd getConcentrations();

        /**
        * @brief Set the alpha coefficients of a 1D-Grid. Grid must be one dimensional. 
        * 
        * @param alpha An Eigen3 MatrixXd with 1 row holding the alpha coefficients. Matrix 
        *              columns must have same size as length of grid. 
        */
        void setAlpha(MatrixXd alpha);

        /**
        * @brief Set the alpha coefficients of a 2D-Grid. Grid must be two dimensional. 
        * 
        * @param alphaX An Eigen3 MatrixXd holding the alpha coefficients in x-direction. 
        *                Matrix must be of same size as the grid. 
        * @param alphaY An Eigen3 MatrixXd holding the alpha coefficients in y-direction.
        *                Matrix must be of same size as the grid. 
        */
        void setAlpha(MatrixXd alphaX, MatrixXd alphaY);

        /**
         * @brief Gets the matrix of alpha coefficients of a 1D-Grid. Grid must be one dimensional.
         * 
         * @return MatrixXd A matrix with 1 row holding the alpha coefficients. 
         */
        const MatrixXd getAlpha();

        /**
         * @brief Gets the matrix of alpha coefficients in x-direction of a 2D-Grid. Grid must be 
         *        two dimensional. 
         * 
         * @return MatrixXd A matrix holding the alpha coefficients in x-direction. 
         */
        const MatrixXd getAlphaX();

        /**
         * @brief Gets the matrix of alpha coefficients in y-direction of a 2D-Grid. Grid must be 
         *        two dimensional. 
         * 
         * @return MatrixXd A matrix holding the alpha coefficients in y-direction. 
         */
        const MatrixXd getAlphaY();

        /**
         * @brief Gets the dimensions of the grid. 
         * 
         * @return int Dimensions, either 1 or 2. 
         */
        int getDim();

        /**
         * @brief Gets length of 1D grid. Must be one dimensional grid.
         * 
         * @return int Length of 1D grid.
         */
        int getLength();

        /**
         * @brief Gets the number of rows of the grid.
         * 
         * @return int Number of rows. 
         */
        int getRow();

        /**
         * @brief Gets the number of columns of the grid.
         * 
         * @return int Number of columns. 
         */
        int getCol();

        /**
         * @brief Sets the domain length of a 1D-Grid. Grid must be one dimensional. 
         * 
         * @param domainLength A double value of the domain length. Must be positive.
         */
        void setDomain(double domainLength);

        /**
         * @brief Sets the domain size of a 2D-Grid. Grid must be two dimensional.
         * 
         * @param domainRow A double value of the domain size in y-direction. Must be positive. 
         * @param domainCol A double value of the domain size in x-direction. Must be positive. 
         */
        void setDomain(double domainRow,double domainCol);

        /**
         * @brief Gets the delta value for 1D-Grid. Grid must be one dimensional.
         * 
         * @return double Delta value. 
         */
        double getDelta();

        /**
         * @brief Gets the delta value in x-direction. 
         * 
         * @return double Delta value in x-direction.
         */
        double getDeltaCol();

        /**
         * @brief Gets the delta value in y-direction. Must be two dimensional grid.
         * 
         * @return double Delta value in y-direction. 
         */
        double getDeltaRow();


    private:

        int col; // number of grid columns
        int row; // number of grid rows
        int dim; // 1D or 2D
        double domainCol; // number of domain columns
        double domainRow; // number of domain rows
        double deltaCol; // delta in x-direction (between columns)
        double deltaRow; // delta in y-direction (between rows)
        MatrixXd concentrations; // Matrix holding grid concentrations 
        MatrixXd alphaX; // Matrix holding alpha coefficients in x-direction
        MatrixXd alphaY; // Matrix holding alpha coefficients in y-direction

};
