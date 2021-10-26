add_library(diffusion OBJECT diffusion.cpp diffusion.hpp)
target_link_libraries(diffusion Eigen3::Eigen)

add_executable(test main.cpp)
target_link_libraries(test PUBLIC diffusion)
