#include <tug/Simulation.hpp>
#include <iostream>
#include <fstream>
#include <chrono>

int main(int argc, char *argv[]) {

    int n[4] = {100, 500, 1000, 2000};
    int threads[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int iterations[1] = {5};
    int repetition = 1;

    ofstream myfile;
    myfile.open("testLarge.csv");

    for (int i = 0; i < size(n); i++){
        cout << "Grid size: " << n[i] << " x " << n[i] << endl << endl;
        myfile << "Grid size: " << n[i] << " x " << n[i] << endl << endl;
        for(int j = 0; j < size(iterations); j++){
            cout << "Iterations: " << iterations[j] << endl;
            myfile << "Iterations: " << iterations[j] << endl;
            for (int k = 0; k < repetition; k++){
                cout << "Wiederholung: " << k << endl;
                Grid grid = Grid(n[i], n[i]);
                grid.setDomain(1, 1);

                MatrixXd concentrations = MatrixXd::Constant(n[i], n[i], 0);
                concentrations(n[i]/2,n[i]/2) = 1;
                grid.setConcentrations(concentrations);
                MatrixXd alpha = MatrixXd::Constant(n[i], n[i], 0.5);

                Boundary bc = Boundary(grid);

                Simulation sim = Simulation(grid, bc, BTCS_APPROACH);
                sim.setSolver(THOMAS_ALGORITHM_SOLVER);

                if(argc == 2){
                    int numThreads = atoi(argv[1]);
                    sim.setNumberThreads(numThreads);
                }
                else{
                    sim.setNumberThreads(1);
                }

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