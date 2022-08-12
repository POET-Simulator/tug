add_library(BTCSBoundary STATIC BTCSBoundaryCondition.cpp)
add_library(BTCSDiffusion STATIC BTCSDiffusion.cpp)

target_include_directories(BTCSBoundary PUBLIC ../include)

target_link_libraries(BTCSDiffusion Eigen3::Eigen BTCSBoundary)

if(BTCS_USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(BTCSDiffusion OpenMP::OpenMP_CXX)
endif()

target_include_directories(BTCSDiffusion PUBLIC ../include)
