#ifndef TUGUTILS_H_
#define TUGUTILS_H_

#include <chrono>
#include <stdexcept>
#include <string>

#define throw_invalid_argument(msg)                                            \
  throw std::invalid_argument(std::string(__FILE__) + ":" +                    \
                              std::to_string(__LINE__) + ":" +                 \
                              std::string(msg))

#define throw_out_of_range(msg)                                                \
  throw std::out_of_range(std::string(__FILE__) + ":" +                        \
                          std::to_string(__LINE__) + ":" + std::string(msg))

#define time_marker() std::chrono::high_resolution_clock::now()

#define diff_time(start, end)                                                  \
  ({                                                                           \
    std::chrono::duration<double> duration =                                   \
        std::chrono::duration_cast<std::chrono::duration<double>>(end -        \
                                                                  start);      \
    duration.count();                                                          \
  })

// calculates arithmetic or harmonic mean of alpha between two cells
constexpr double calcAlphaIntercell(double alpha1, double alpha2,
                                    bool useHarmonic = true) {
  if (useHarmonic) {
    return double(2) / ((double(1) / alpha1) + (double(1) / alpha2));
  } else {
    return 0.5 * (alpha1 + alpha2);
  }
}
#endif // TUGUTILS_H_
