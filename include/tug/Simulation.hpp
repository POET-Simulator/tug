/**
 * @file Simulation.hpp
 * @brief API of Simulation class, that holds all information regarding a specific simulation
 *        run like its timestep, number of iterations and output options. Simulation object
 *        also holds a predefined Grid and Boundary object. 
 *
 */
#include "Boundary.hpp"
#include <ios>

using namespace std;

/**
 * @brief Enum defining the two implemented solution approaches. 
 * 
 */
enum APPROACH {
    FTCS_APPROACH, // Forward Time-Centered Space
    BTCS_APPROACH, // Backward Time-Centered Space solved with EigenLU solver
    CRANK_NICOLSON_APPROACH 
};

/**
 * @brief Enum defining the Linear Equation solvers
 * 
 */
enum SOLVER {
    EIGEN_LU_SOLVER, // EigenLU solver
    THOMAS_ALGORITHM_SOLVER // Thomas Algorithm solver; more efficient for tridiagonal matrices
};

/**
 * @brief Enum holding different options for .csv output.
 * 
 */
enum CSV_OUTPUT {
    CSV_OUTPUT_OFF, // do not produce csv output
    CSV_OUTPUT_ON, // produce csv output with last concentration matrix
    CSV_OUTPUT_VERBOSE, // produce csv output with all concentration matrices
    CSV_OUTPUT_XTREME // csv output like VERBOSE but additional boundary conditions at beginning
};

/**
 * @brief Enum holding different options for console output.
 * 
 */
enum CONSOLE_OUTPUT {
    CONSOLE_OUTPUT_OFF, // do not print any output to console
    CONSOLE_OUTPUT_ON, // print before and after concentrations to console
    CONSOLE_OUTPUT_VERBOSE // print all concentration matrices to console
};

/**
 * @brief Enum holding different options for time measurement. 
 * 
 */
enum TIME_MEASURE {
    TIME_MEASURE_OFF, // do not print any time measures
    TIME_MEASURE_ON // print time measure after last iteration
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
       * @brief Set up a simulation environment. The timestep and number of iterations
       *        must be set. For the BTCS approach, the Thomas algorithm is used as 
       *        the default linear equation solver as this is faster for tridiagonal
       *        matrices. CSV output, console output and time measure are off by default. 
       *        Also, the number of cores is set to the maximum number of cores by default.
       *
       * @param grid Valid grid object
       * @param bc Valid boundary condition object
       * @param approach Approach to solving the problem. Either FTCS or BTCS.
       */
      Simulation(Grid &grid, Boundary &bc, APPROACH approach);

      /**
       * @brief Set the option to output the results to a CSV file. Off by default.
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
       * @brief Set the options for outputting information to the console. Off by default. 
       *
       * @param console_output Valid output option. The following options can be set
       *                       here:
       *                        - CONSOLE_OUTPUT_OFF: do not print any output to console
       *                        - CONSOLE_OUTPUT_ON: print before and after concentrations to console
       *                        - CONSOLE_OUTPUT_VERBOSE: print all concentration matrices to console
       */
      void setOutputConsole(CONSOLE_OUTPUT console_output);

      /**
       * @brief Set the Time Measure option. Off by default. 
       *
       * @param time_measure The following options are allowed:
       *                     - TIME_MEASURE_OFF: Time of simulation is not printed to console
       *                     - TIME_MEASURE_ON: Time of simulation run is printed to console
       */
      void setTimeMeasure(TIME_MEASURE time_measure);

      /**
       * @brief Setting the time step for each iteration step. Time step must be
       *        greater than zero. Setting the timestep is required. 
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
       *        than zero must be specified here. Setting iterations is required. 
       *
       * @param iterations Number of iterations to be simulated.
       */
      void setIterations(int iterations);

      /**
       * @brief Set the desired linear equation solver to be used for BTCS approach. Without effect
       *        in case of FTCS approach.
       * 
       * @param solver Solver to be used. Default is Thomas Algorithm as it is more efficient for 
       *               tridiagonal Matrices. 
       */
      void setSolver(SOLVER solver);

      /**
       * @brief Set the number of desired openMP Threads.
       *
       * @param num_threads Number of desired threads. Must have a value between
       *                    1 and the maximum available number of processors. The maximum number of
       *                    processors is set as the default case during Simulation construction.
       */
      void setNumberThreads(int num_threads);

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
       *        <approach> + <number rows> + <number columns> + <number of iterations>+<counter>.csv 
       *
       * @return string Filename with configured simulation parameters.
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
        int innerIterations;
        int numThreads;
        CSV_OUTPUT csv_output;
        CONSOLE_OUTPUT console_output;
        TIME_MEASURE time_measure;

        Grid &grid;
        Boundary &bc;
        APPROACH approach;
        SOLVER solver;

};
