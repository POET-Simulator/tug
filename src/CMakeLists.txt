add_library(diffusion OBJECT diffusion.cpp diffusion.hpp)
target_link_libraries(diffusion Eigen3::Eigen)

add_library(diffusion_class OBJECT BTCSDiffusion.cpp BTCSDiffusion.hpp)
target_link_libraries(diffusion_class Eigen3::Eigen)

add_executable(test main.cpp)
target_link_libraries(test PUBLIC diffusion)

add_executable(test_class test_class.cpp)
target_link_libraries(test_class PUBLIC diffusion_class)
