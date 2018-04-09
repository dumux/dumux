add_subdirectory("thermal")
#install headers
install(FILES
localresidual.hh
gridvariables.hh
indices.hh
volumevariables.hh
vtkoutputfields.hh
model.hh
newtonsolver.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/nonequilibrium)
