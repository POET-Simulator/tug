#include <Eigen/Core>

using namespace Eigen;

class Grid {
    public:

        /**
        * @brief Construct a new Grid object
        * 
        * @param n 
        */
        Grid(int n);

        /**
        * @brief Construct a new Grid object
        * 
        * @param n 
        * @param m 
        */
        Grid(int n, int m);

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


    private:

        int dim;
        int n;
        int m;
        Matrix2d concentrations;
        Matrix2d alpha_x;
        Matrix2d alpha_y;
};