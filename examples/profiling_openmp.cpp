#include <tug/Simulation.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <easy/profiler.h>


int main(int argc, char *argv[]) {
    EASY_MAIN_THREAD;
    EASY_PROFILER_ENABLE;
    profiler::startListen();  

    int n = 1000;

    Grid grid = Grid(n, n);
    grid.setDomain(10, 10);

    MatrixXd concentrations = MatrixXd::Constant(n, n, 0);
    concentrations(n/2,n/2) = 1;
    grid.setConcentrations(concentrations);
    MatrixXd alpha = MatrixXd::Constant(n, n, 0.001);
    
    Boundary bc = Boundary(grid);
    
    Simulation sim = Simulation(grid, bc, BTCS_APPROACH);
    sim.setSolver(THOMAS_ALGORITHM_SOLVER);
    sim.setNumberThreads(1);
    
    sim.setTimestep(0.001);
    sim.setIterations(2);
    sim.setOutputCSV(CSV_OUTPUT_OFF);
    
    auto begin = std::chrono::high_resolution_clock::now();

    EASY_BLOCK("SIMULATION");
    sim.run();
    EASY_END_BLOCK;

    auto end = std::chrono::high_resolution_clock::now();

    auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    cout << milliseconds.count() << endl;

    profiler::dumpBlocksToFile("./mytest_profile.prof");
    profiler::stopListen();

    return(0);
}
