#pragma once

#include <cassert>

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

#define tug_assert(expr, msg) assert((expr) && msg)

// calculates arithmetic or harmonic mean of alpha between two cells
template <typename T>
constexpr T calcAlphaIntercell(T alpha1, T alpha2, bool useHarmonic = true) {
  if (useHarmonic) {
    const T operand1 = alpha1 == 0 ? 0 : 1 / alpha1;
    const T operand2 = alpha2 == 0 ? 0 : 1 / alpha2;

    const T denom = operand1 + operand2;

    return denom == 0 ? 0 : 2. / denom;
  } else {
    return 0.5 * (alpha1 + alpha2);
  }
}
