add_executable(naaice BTCS_2D_NAAICE.cpp)

get_filename_component(IN_CONC_FILE "init_conc.csv" REALPATH)
get_filename_component(IN_ALPHAX_FILE "alphax.csv" REALPATH)

configure_file(files.hpp.in files.hpp)
target_include_directories(naaice PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(naaice PUBLIC tug)
