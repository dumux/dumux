add_input_file_links()

add_dumux_test(tutorial_sequential tutorial_sequential tutorial_sequential.cc
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial_sequential)

add_dumux_test(tutorial_implicit tutorial_implicit tutorial_implicit.cc
  ${CMAKE_CURRENT_BINARY_DIR}/tutorial_implicit)

#install sources
install(FILES
tutorial_implicit.cc
tutorial_sequential.cc
tutorialproblem_implicit.hh
tutorialproblem_sequential.hh
tutorialspatialparams_implicit.hh
tutorialspatialparams_sequential.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/tutorial)
add_subdirectory(dumux-course)
