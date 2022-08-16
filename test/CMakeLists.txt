include(FetchContent)

FetchContent_Declare(
    DocTest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.9
)

FetchContent_MakeAvailable(DocTest)

#add_library(doctest INTERFACE)
#target_include_directories(doctest INTERFACE doctest)

add_executable(test setup.cpp testBoundaryCondition.cpp testDiffusion.cpp)
target_link_libraries(test doctest BTCSDiffusion)
