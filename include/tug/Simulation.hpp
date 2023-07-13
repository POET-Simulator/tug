#include "Boundary.hpp"

enum APPROACH {
    FTCS,
    BTCS
};

enum CSV_OUTPUT {
    CSV_OUTPUT_OFF,
    CSV_OUTPUT_ON,
    CSV_OUTPUT_VERBOSE
};

class Simulation {
    public:

        /**
         * @brief Construct a new Simulation object
         * 
         * @param grid 
         * @param bc 
         * @param aproach 
         */
        Simulation(Grid grid, Boundary bc, APPROACH aproach);

        /**
         * @brief 
         * 
         * @param csv_output 
         */
        void outputCSV(CSV_OUTPUT csv_output);

        /**
         * @brief Set the Timestep object
         * 
         * @param timetstep 
         */
        void setTimestep(double timetstep);

        /**
         * @brief Get the Timestep object
         * 
         */
        void getTimestep();

        /**
         * @brief Set the Iterations object
         * 
         * @param iterations 
         */
        void setIterations(int iterations);

        /**
         * @brief Get the Iterations object
         * 
         * @return auto 
         */
        auto getIterations();

        /**
         * @brief 
         * 
         * @return auto 
         */
        auto run();

    private:

        double timestep;
        int iterations;
        CSV_OUTPUT csv_output;

        Boundary bc;
        APPROACH approach;

};
