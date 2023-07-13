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
        auto getConcentrations();

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

        int getDim();

        auto getRow();

        auto getCol();


    private:

        int dim;
        int row;
        int col;
        Matrix2d concentrations;
        Matrix2d alpha_x;
        Matrix2d alpha_y;
};