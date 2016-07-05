add_subdirectory("tpfa")

install(FILES
fluxvariablescachevector.hh
stencils.hh
volumevariablesvector.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered)