add_input_file_links()

add_dumux_test(tutorial_implicit tutorial_implicit tutorial_implicit.cc
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial_implicit)

#install sources
install(FILES
tutorial_implicit.cc
tutorialproblem_implicit.hh
tutorialspatialparams_implicit.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/tutorial/tutorial_implicit)
