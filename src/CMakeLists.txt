add_library(tug Boundary.cpp Grid.cpp Simulation.cpp FTCS.cpp BTCSv2.cpp)

target_link_libraries(tug Eigen3::Eigen)

if(TUG_USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(tug OpenMP::OpenMP_CXX)
endif()

target_include_directories(tug PUBLIC ../include)
