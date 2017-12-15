add_subdirectory("sequential")

#install headers
install(FILES
indices.hh
incompressiblelocalresidual.hh
model.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/1p)
