add_subdirectory("cellcentered")

#install headers
install(FILES
problem.hh
properties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/1p/sequential/diffusion)
