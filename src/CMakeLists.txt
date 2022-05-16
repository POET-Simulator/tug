set(HEADER_LIST "${BTCSDiffusion_SOURCE_DIR}/include/diffusion/BTCSDiffusion.hpp"
  "${BTCSDiffusion_SOURCE_DIR}/include/diffusion/BoundaryCondition.hpp")

add_library(BTCSDiffusion STATIC BTCSDiffusion.cpp ${HEADER_LIST})

target_link_libraries(BTCSDiffusion Eigen3::Eigen)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(BTCSDiffusion OpenMP::OpenMP_CXX)
endif()

target_include_directories(BTCSDiffusion PUBLIC ../include)
