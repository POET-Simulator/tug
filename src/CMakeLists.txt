add_library(tug BTCSDiffusion.cpp grid/BoundaryCondition.cpp)

target_link_libraries(tug Eigen3::Eigen)

if(TUG_USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(tug OpenMP::OpenMP_CXX)
endif()

target_include_directories(tug PUBLIC ../include)
