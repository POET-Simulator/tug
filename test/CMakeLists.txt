add_library(doctest INTERFACE)
target_include_directories(doctest INTERFACE doctest)

add_executable(test setup.cpp testBoundaryCondition.cpp testDiffusion.cpp)
target_link_libraries(test doctest BTCSDiffusion)
