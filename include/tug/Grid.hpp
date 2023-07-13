#include <iostream>
#include <vector>

using namespace std;

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
     * @brief Set the Alpha object
     * 
     * @param alpha 
     */
    void setAlpha(vector<float> alpha);

    void setAlpha(vector<float> alpha_x, vector<float> alpha_y);

    private:
    

};