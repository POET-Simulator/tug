add_library(diffusion OBJECT BTCSDiffusion.cpp BTCSDiffusion.hpp)
target_link_libraries(diffusion Eigen3::Eigen)

add_executable(1D main_1D.cpp)
target_link_libraries(1D PUBLIC diffusion)

add_executable(2D main_2D.cpp)
target_link_libraries(2D PUBLIC diffusion)

add_executable(Comp2D main_2D_mdl.cpp)
target_link_libraries(Comp2D PUBLIC diffusion)
