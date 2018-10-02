add_subdirectory(staggered)

install(FILES
fluxvariables.hh
fluxvariablescache.hh
indices.hh
localresidual.hh
model.hh
problem.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/navierstokes)
