find_package(OpenMP)

add_library(diffusion OBJECT BTCSDiffusion.cpp BTCSDiffusion.hpp)
target_link_libraries(diffusion Eigen3::Eigen OpenMP::OpenMP_CXX)
target_include_directories(diffusion PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
