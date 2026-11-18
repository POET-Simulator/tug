# add_executable(BTCS_2D_proto_example BTCS_2D_proto_example.cpp)
# add_executable(FTCS_2D_proto_example_mdl FTCS_2D_proto_example_mdl.cpp)
# add_executable(FTCS_2D_proto_closed_mdl FTCS_2D_proto_closed_mdl.cpp)
add_executable(Advection Advection.cpp)

# target_link_libraries(BTCS_2D_proto_example tug)
# target_link_libraries(FTCS_2D_proto_closed_mdl tug)
# target_link_libraries(FTCS_2D_proto_example_mdl tug)
target_link_libraries(Advection tug)
