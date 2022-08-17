add_library(BTCSDiffusion BTCSDiffusion.cpp grid/BoundaryCondition.cpp)

target_link_libraries(BTCSDiffusion Eigen3::Eigen)

if(BTCS_USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(BTCSDiffusion OpenMP::OpenMP_CXX)
endif()

target_include_directories(BTCSDiffusion PUBLIC ../include)
