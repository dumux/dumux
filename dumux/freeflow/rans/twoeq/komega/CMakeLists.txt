add_subdirectory(staggered)

install(FILES
fluxvariables.hh
indices.hh
localresidual.hh
model.hh
models.hh
problem.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans/twoeq/komega)
