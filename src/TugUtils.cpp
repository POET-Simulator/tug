#include <chrono>
#include <fstream>
#include <stdexcept>
#include <string>

using namespace std;

// used for throwing an invalid argument message
#define throw_invalid_argument(msg)                                            \
  throw std::invalid_argument(std::string(__FILE__) + ":" +                    \
                              std::to_string(__LINE__) + ":" +                 \
                              std::string(msg))

// used for throwing an out of range message
#define throw_out_of_range(msg)                                                \
  throw std::out_of_range(std::string(__FILE__) + ":" +                        \
                          std::to_string(__LINE__) + ":" + std::string(msg))

// get current time
#define time_marker() std::chrono::high_resolution_clock::now()

// calculates difference between two time points
#define diff_time(start, end)                                                  \
  ({                                                                           \
    std::chrono::duration<double> duration =                                   \
        std::chrono::duration_cast<std::chrono::duration<double>>(end -        \
                                                                  start);      \
    duration.count();                                                          \
  })

// calculates arithmetic or harmonic mean of alpha between two cells
static double calcAlphaIntercell(const double &alpha1, const double &alpha2,
                                 bool useHarmonic = true) {
  if (useHarmonic) {
    return double(2) / ((double(1) / alpha1) + (double(1) / alpha2));
  } else {
    return 0.5 * (alpha1 + alpha2);
  }
}
