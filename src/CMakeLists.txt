set(HEADER_LIST "${Diffusion_SOURCE_DIR}/include/diffusion/BTCSDiffusion.hpp"
  "${Diffusion_SOURCE_DIR}/include/diffusion/BoundaryCondition.hpp")

add_library(diffusion STATIC BTCSDiffusion.cpp ${HEADER_LIST})

target_link_libraries(diffusion Eigen3::Eigen)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(diffusion OpenMP::OpenMP_CXX)
endif()

target_include_directories(diffusion PUBLIC ../include)
