add_subdirectory("sequential")

#install headers
install(FILES
adaptionhelper.hh
gridadaptindicator.hh
indices.hh
incompressiblelocalresidual.hh
model.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p)
