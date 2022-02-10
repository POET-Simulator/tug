add_library(diffusion OBJECT BTCSDiffusion.cpp BTCSDiffusion.hpp)
target_link_libraries(diffusion Eigen3::Eigen)

add_executable(test main.cpp)
target_link_libraries(test PUBLIC diffusion)

add_executable(2D main_2D.cpp)
target_link_libraries(2D PUBLIC diffusion)
