add_library(diffusion OBJECT BTCSDiffusion.cpp BTCSDiffusion.hpp)
target_link_libraries(diffusion Eigen3::Eigen)
target_include_directories(diffusion PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
