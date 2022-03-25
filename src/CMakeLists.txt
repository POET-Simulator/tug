find_package(OpenMP)

add_library(diffusion OBJECT BTCSDiffusion.cpp BTCSDiffusion.hpp)

target_link_libraries(diffusion Eigen3::Eigen)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(diffusion OpenMP::OpenMP_CXX)
endif()

target_include_directories(diffusion PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
