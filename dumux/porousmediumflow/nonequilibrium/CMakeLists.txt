add_subdirectory("thermal")
#install headers
install(FILES
localresidual.hh
gridvariables.hh
indices.hh
newtoncontroller.hh
volumevariables.hh
vtkoutputfields.hh
model.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/nonequilibrium)
