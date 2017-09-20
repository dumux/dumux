add_input_file_links()

add_dumux_test(ex3_a ex3_a ex3_a.cc  ${CMAKE_CURRENT_BINARY_DIR}/exercise3)
add_dumux_test(ex3_b ex3_b ex3_b.cc  ${CMAKE_CURRENT_BINARY_DIR}/exercise3)

#install sources
install(FILES
exercise3.cc
exercise3_problem.hh
exercise3_spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/tutorial/material_exercise)
