add_subdirectory(sequential)

install(FILES
indices.hh
model.hh
primaryvariableswitch.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p2c)
