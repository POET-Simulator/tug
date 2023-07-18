#include <Eigen/Core>

using namespace Eigen;

class Grid {
    public:

        /**
        * @brief Construct a new Grid object
        * 
        * @param col 
        */
        Grid(int col);

        /**
        * @brief Construct a new Grid object
        * 
        * @param row 
        * @param col 
        */
        Grid(int row, int col);

        /**
        * @brief Set the Concentrations object
        * 
        * @param concentrations 
        */
        void setConcentrations(Matrix2d concentrations);

        /**
         * @brief Get the Concentrations object
         * 
         * @return auto 
         */
        Matrix2d getConcentrations();

        /**
        * @brief Set the Alpha object
        * 
        * @param alpha 
        */
        void setAlpha(Matrix2d alpha);

        /**
        * @brief Set the Alpha object
        * 
        * @param alpha_x 
        * @param alpha_y 
        */
        void setAlpha(Matrix2d alpha_x, Matrix2d alpha_y);

        Matrix2d getAlphaX();

        Matrix2d getAlphaY();

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
        Matrix2d concentrations;
        Matrix2d alpha_x;
        Matrix2d alpha_y;
};