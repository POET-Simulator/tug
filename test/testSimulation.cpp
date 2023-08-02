#include <stdio.h>
#include <doctest/doctest.h>
#include <tug/Simulation.hpp>
#include <filesystem>
#include "TestUtils.cpp"
#include <string>

namespace fs = std::filesystem;

static Grid setupSimulation() {
    int row = 11;
    int col = 11;
    int domain_row = 10;
    int domain_col = 10;


    // Grid
    Grid grid = Grid(row, col);
    grid.setDomain(domain_row, domain_col);

    MatrixXd concentrations = MatrixXd::Constant(row, col, 0);
    concentrations(5,5) = 1;
    grid.setConcentrations(concentrations);

    MatrixXd alpha = MatrixXd::Constant(row, col, 1);
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 6; j++) {
            alpha(i, j) = 0.01;
        }
    }
    for (int i = 0; i < 5; i++) {
        for (int j = 6; j < 11; j++) {
            alpha(i, j) = 0.001;
        }
    }
    for (int i = 5; i < 11; i++) {
        for (int j = 6; j < 11; j++) {
            alpha(i, j) = 0.1;
        }
    }
    grid.setAlpha(alpha, alpha);


    // Boundary
    Boundary bc = Boundary(grid);


    // Simulation
    Simulation sim = Simulation(grid, bc, FTCS_APPROACH);
    sim.setTimestep(0.001);
    sim.setIterations(7000);
    // sim.setOutputCSV(CSV_OUTPUT_ON);
    // sim.setOutputConsole(CONSOLE_OUTPUT_ON);


    // RUN
    return sim.run();
}

TEST_CASE("equality to reference matrix") {
    fs::path p = fs::current_path().parent_path().parent_path();
    string test_path = p.generic_string() + "/test/FTCS_11_11_7000.csv";
    MatrixXd reference = CSV2Eigen(test_path);

    Grid grid = setupSimulation();
    CHECK(checkSimilarity(reference, grid.getConcentrations(), 0.1) == true);
}