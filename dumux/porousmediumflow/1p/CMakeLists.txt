add_subdirectory(sequential)

install(FILES
incompressiblelocalresidual.hh
indices.hh
model.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/1p)
