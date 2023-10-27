#include <Eigen/Eigen>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <tug/Simulation.hpp>

using namespace Eigen;
using namespace std;
using namespace tug;

int main(int argc, char *argv[]) {

  int n[] = {2000};
  int threads[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  int iterations[1] = {1};
  int repetition = 10;

  for (int l = 0; l < size(threads); l++) {
    // string filename = "ftcs_openmp_" + to_string(threads[l]) + ".csv";
    ofstream myfile;
    myfile.open("speedup_1000.csv", std::ios::app);
    myfile << "Number threads: " << threads[l] << endl;

    for (int i = 0; i < size(n); i++) {
      cout << "Grid size: " << n[i] << " x " << n[i] << endl << endl;
      // myfile << "Grid size: " << n[i] << " x " << n[i] << endl << endl;
      for (int j = 0; j < size(iterations); j++) {
        cout << "Iterations: " << iterations[j] << endl;
        // myfile << "Iterations: " << iterations[j] << endl;
        for (int k = 0; k < repetition; k++) {
          cout << "Wiederholung: " << k << endl;
          Grid64 grid(n[i], n[i]);
          grid.setDomain(1, 1);

          MatrixXd concentrations = MatrixXd::Constant(n[i], n[i], 0);
          concentrations(n[i] / 2, n[i] / 2) = 1;
          grid.setConcentrations(concentrations);
          MatrixXd alpha = MatrixXd::Constant(n[i], n[i], 0.5);

          Boundary bc = Boundary(grid);

          Simulation sim = Simulation(grid, bc);

          if (argc == 2) {
            int numThreads = atoi(argv[1]);
            sim.setNumberThreads(numThreads);
          } else {
            sim.setNumberThreads(threads[l]);
          }

          sim.setTimestep(0.01);
          sim.setIterations(iterations[j]);
          sim.setOutputCSV(CSV_OUTPUT_OFF);

          auto begin = std::chrono::high_resolution_clock::now();
          sim.run();
          auto end = std::chrono::high_resolution_clock::now();
          auto milliseconds =
              std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                    begin);
          myfile << milliseconds.count() << endl;
        }
      }
      cout << endl;
      myfile << endl;
    }
    myfile.close();
  }
}
