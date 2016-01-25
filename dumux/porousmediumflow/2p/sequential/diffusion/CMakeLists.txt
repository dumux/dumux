add_subdirectory("cellcentered")
add_subdirectory("mpfa")
add_subdirectory("mimetic")

#install headers
install(FILES
problem.hh
properties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential/diffusion)
