add_subdirectory("staggered")

#install headers
install(FILES
fluxvariables.hh
indices.hh
localresidual.hh
navierstokesncmodel.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/compositional)
