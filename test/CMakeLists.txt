include(FetchContent)

FetchContent_Declare(
    DocTest
    GIT_REPOSITORY https://github.com/doctest/doctest.git
    GIT_TAG v2.4.9
)

FetchContent_MakeAvailable(DocTest)

add_executable(testTug setup.cpp testBoundaryCondition.cpp testDiffusion.cpp)
target_link_libraries(testTug doctest tug)

add_custom_target(check
    COMMAND $<TARGET_FILE:testTug>
    DEPENDS testTug
)
