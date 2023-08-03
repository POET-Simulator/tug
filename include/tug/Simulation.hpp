#include "Boundary.hpp"
#include <ios>

using namespace std;

enum APPROACH {
    FTCS_APPROACH, // Forward Time-Centered Space
    BTCS_APPROACH // Backward Time-Centered Space
};

enum CSV_OUTPUT {
    CSV_OUTPUT_OFF, // do not produce csv output
    CSV_OUTPUT_ON, // produce csv output with last concentration matrix
    CSV_OUTPUT_VERBOSE, // produce csv output with all concentration matrices
    CSV_OUTPUT_XTREME // produce csv output with all concentration matrices and simulation environment
};

enum CONSOLE_OUTPUT {
    CONSOLE_OUTPUT_OFF, // do not print any output to console
    CONSOLE_OUTPUT_ON, // print before and after concentrations to console
    CONSOLE_OUTPUT_VERBOSE // print all concentration matrices to console
};

enum TIME_MEASURE {
    TIME_MEASURE_OFF, // do not print any time measures
    TIME_MEASURE_ON, // print time measure after last iteration
    TIME_MEASURE_VERBOSE // print time measures after each iteration
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
        Simulation(Grid &grid, Boundary &bc, APPROACH approach);

        /**
         * @brief 
         * 
         * @param csv_output 
         */
        void setOutputCSV(CSV_OUTPUT csv_output);

        /**
         * @brief Set the Output Console object
         * 
         * @param console_output 
         */
        void setOutputConsole(CONSOLE_OUTPUT console_output);

        /**
         * @brief Set the Time Measure object
         * 
         * @param time_measure 
         */
        void setTimeMeasure(TIME_MEASURE time_measure);

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
        double getTimestep();

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
        int getIterations();

        /**
         * @brief Print the current concentrations of the grid to standard out. 
         * 
         */
        void printConcentrationsConsole();

        /**
         * @brief 
         * 
         * @return string 
         */
        string createCSVfile();

        void printConcentrationsCSV(string filename);

        /**
         * @brief 
         * 
         * @return Grid 
         */
        void run();

    private:

        double timestep;
        int iterations;
        CSV_OUTPUT csv_output;
        CONSOLE_OUTPUT console_output;
        TIME_MEASURE time_measure;

        Grid &grid;
        Boundary &bc;
        APPROACH approach;

};
