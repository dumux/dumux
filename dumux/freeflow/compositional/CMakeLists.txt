add_subdirectory("staggered")

#install headers
install(FILES
fluxvariables.hh
indices.hh
localresidual.hh
navierstokesncmodel.hh
lowrekepsilonncmodel.hh
volumevariables.hh
vtkoutputfields.hh
zeroeqncmodel.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/compositional)
