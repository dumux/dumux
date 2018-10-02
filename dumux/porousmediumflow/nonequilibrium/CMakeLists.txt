add_subdirectory(thermal)

install(FILES
gridvariables.hh
indices.hh
localresidual.hh
model.hh
newtonsolver.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/nonequilibrium)
