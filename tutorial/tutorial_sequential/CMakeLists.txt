add_input_file_links()

add_dumux_test(tutorial_sequential tutorial_sequential tutorial_sequential.cc
              ${CMAKE_CURRENT_BINARY_DIR}/tutorial_sequential)

#install sources
install(FILES
tutorial_sequential.cc
tutorialproblem_sequential.hh
tutorialspatialparams_sequential.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/tutorial/tutorial_sequential)
