#include "Boundary.hpp"

using namespace std;

enum APPROACH {
    FTCS_APPROACH,
    BTCS_APPROACH
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
        Simulation(Grid grid, Boundary bc, APPROACH approach);

        /**
         * @brief 
         * 
         * @param csv_output 
         */
        void setOutputCSV(CSV_OUTPUT csv_output);

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
        auto getTimestep();

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

        void printConcentrationsConsole();

        void printConcentrationsCSV(string ident);

        /**
         * @brief 
         * 
         * @return auto 
         */
        void run();

    private:

        double timestep;
        int iterations;
        CSV_OUTPUT csv_output;

        Grid grid;
        Boundary bc;
        APPROACH approach;

};
