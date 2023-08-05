/**
 * @file Simulation.hpp
 * @brief  
 */
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

/**
 * @brief The class forms the interface for performing the diffusion simulations
 * and contains all the methods for controlling the desired parameters, such as
 * time step, number of simulations, etc.
 *
 */
class Simulation {
    public:
      /**
       * @brief Set up a runnable simulation environment with the largest stable
       *        time step and 1000 iterations by passing the required parameters.
       *
       * @param grid Valid grid object
       * @param bc Valid boundary condition object
       * @param approach Approach to solving the problem. Either FTCS or BTCS.
       */
      Simulation(Grid &grid, Boundary &bc, APPROACH approach);

      /**
       * @brief Set the option to output the results to a CSV file.
       *
       *
       * @param csv_output Valid output option. The following options can be set
       *                   here:
       *                     - CSV_OUTPUT_OFF: do not produce csv output
       *                     - CSV_OUTPUT_ON: produce csv output with last
       *                       concentration matrix
       *                     - CSV_OUTPUT_VERBOSE: produce csv output with all
       *                       concentration matrices
       *                     - CSV_OUTPUT_XTREME: produce csv output with all
       *                       concentration matrices and simulation environment
       */
      void setOutputCSV(CSV_OUTPUT csv_output);

      /**
       * @brief Set the options for outputting information to the console.
       *
       * @param console_output Valid output option. The following options can be set
       *                       here:
       *                        - CONSOLE_OUTPUT_OFF: do not print any output to console
       *                        - CONSOLE_OUTPUT_ON: print before and after concentrations to console
       *                        - CONSOLE_OUTPUT_VERBOSE: print all concentration matrices to console
       */
      void setOutputConsole(CONSOLE_OUTPUT console_output);

      /**
       * @brief Set the Time Measure object
       *
       * @param time_measure
       */
      void setTimeMeasure(TIME_MEASURE time_measure);

      /**
       * @brief Setting the time step for each iteration step. Time step must be
       *        greater than zero.
       *
       * @param timestep Valid timestep greater than zero. 
       */
      void setTimestep(double timestep);

      /**
       * @brief Currently set time step is returned.
       * 
       * @return double timestep
       */
      double getTimestep();

      /**
       * @brief Set the desired iterations to be calculated. A value greater
       *        than zero must be specified here.
       *
       * @param iterations Number of iterations to be simulated.
       */
      void setIterations(int iterations);

      /**
       * @brief Return the currently set iterations to be calculated.
       * 
       * @return int Number of iterations.
       */
      int getIterations();

      /**
       * @brief Outputs the current concentrations of the grid on the console.
       *
       */
      void printConcentrationsConsole();

      /**
       * @brief Creates a CSV file with a name containing the current simulation
       *        parameters. If the data name already exists, an additional counter is
       *        appended to the name. The name of the file is built up as follows:
       *        <approach> + <number rows> + <number columns> + <number of iterations>-<counter>.csv 
       *
       * @return string Filename with given simulation parameter.
       */
      string createCSVfile();

      /**
       * @brief Writes the currently calculated concentration values of the grid
       *        into the CSV file with the passed filename.
       *
       * @param filename Name of the file to which the concentration values are
       *                 to be written.
       */
      void printConcentrationsCSV(string filename);

      /**
       * @brief Method starts the simulation process with the previously set
       *        parameters.
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
