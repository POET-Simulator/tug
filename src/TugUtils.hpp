#ifndef BTCSUTILS_H_
#define BTCSUTILS_H_

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
#endif // BTCSUTILS_H_
