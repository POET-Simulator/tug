#include <Eigen/Core>

using namespace Eigen;

class Grid {
    public:

        /**
         * @brief Construct a new 1D Grid object, which represents the fields holding the concentrations.
         * 
         * @param col Defines the length of the grid in terms of cells. Input must be an integer > 3. 
         */
        Grid(int col);

        /**
        * @brief Construct a new Grid object
        * 
        * @param row 
        * @param col 
        */
        Grid(int row, int col);

        // TODO
        // Grid(MatrixXd concentrations);

        /**
        * @brief Set the Concentrations object
        * 
        * @param concentrations 
        */
        void setConcentrations(MatrixXd concentrations);

        /**
         * @brief Get the Concentrations object
         * 
         * @return MatrixXd 
         */
        MatrixXd getConcentrations();

        /**
        * @brief Set the Alpha object
        * 
        * @param alpha 
        */
        void setAlpha(MatrixXd alpha);

        /**
        * @brief Set the Alpha object
        * 
        * @param alpha_x 
        * @param alpha_y 
        */
        void setAlpha(MatrixXd alpha_x, MatrixXd alpha_y);

        MatrixXd getAlphaX();

        MatrixXd getAlphaY();

        int getDim();

        int getRow();

        int getCol();

        void setDomain(int domain_col);

        void setDomain(int domain_row, int domain_col);

        double getDeltaCol();

        double getDeltaRow();


    private:

        int dim;
        int col;
        int row;
        int domain_col;
        int domain_row;
        double delta_col;
        double delta_row;
        MatrixXd concentrations;
        MatrixXd alpha_x;
        MatrixXd alpha_y;
};