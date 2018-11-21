add_subdirectory(thermal)

install(FILES
gridvariables.hh
indices.hh
iofields.hh
localresidual.hh
model.hh
newtonsolver.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/nonequilibrium)
