#include <tug/Simulation.hpp>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {

    int n[4] = {10, 20, 50, 100};
    int iterations[1] = {100};
    int repetition = 10;

    ofstream myfile;
    myfile.open("btcs_time_measurement_openmp_10.csv");

    for (int i = 0; i < size(n); i++){
        cout << "Grid size: " << n[i] << " x " << n[i] << endl << endl;
        myfile << "Grid size: " << n[i] << " x " << n[i] << endl;
        for(int j = 0; j < size(iterations); j++){
            cout << "Iterations: " << iterations[j] << endl;
            for (int k = 0; k < repetition; k++){
                cout << "Wiederholung: " << k << endl;
                Grid grid = Grid(n[i], n[i]);
                grid.setDomain(n[i], n[i]);

                MatrixXd concentrations = MatrixXd::Constant(n[i], n[i], 0);
                concentrations(5,5) = 1;
                grid.setConcentrations(concentrations);
                MatrixXd alpha = MatrixXd::Constant(n[i], n[i], 0.5);

                Boundary bc = Boundary(grid);

                Simulation sim = Simulation(grid, bc, BTCS_APPROACH);


                sim.setTimestep(0.001);
                sim.setIterations(iterations[j]);
                sim.setOutputCSV(CSV_OUTPUT_OFF);

                auto begin = std::chrono::high_resolution_clock::now();
                sim.run();
                auto end = std::chrono::high_resolution_clock::now();
                auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
                myfile << milliseconds.count() << endl;
            }
        }
        cout << endl;
        myfile << endl;

    }
    myfile.close();
}