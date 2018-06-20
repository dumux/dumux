add_subdirectory(cellcentered)
add_subdirectory(mimetic)
add_subdirectory(mpfa)

install(FILES
problem.hh
properties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential/diffusion)
