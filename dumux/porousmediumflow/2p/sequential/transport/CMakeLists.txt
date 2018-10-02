add_subdirectory(cellcentered)

install(FILES
problem.hh
properties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential/transport)
