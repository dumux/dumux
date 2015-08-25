add_input_file_links()

add_dumux_test(tutorial_decoupled tutorial_decoupled tutorial_decoupled.cc
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial_decoupled)

add_dumux_test(tutorial_coupled tutorial_coupled tutorial_coupled.cc
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial_coupled)

#install sources
install(FILES
tutorial_coupled.cc
tutorial_decoupled.cc
tutorialproblem_coupled.hh
tutorialproblem_decoupled.hh
tutorialspatialparams_coupled.hh
tutorialspatialparams_decoupled.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/tutorial)
