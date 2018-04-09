
#install headers
install(FILES
indices.hh
localresidual.hh
model.hh
newtonsolver.hh
privarswitchnewtonsolver.hh
primaryvariableswitch.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/richards)
