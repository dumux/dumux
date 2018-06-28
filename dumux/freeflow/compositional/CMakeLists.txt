add_subdirectory(staggered)

install(FILES
fluxvariables.hh
indices.hh
kepsilonncmodel.hh
komegancmodel.hh
localresidual.hh
lowrekepsilonncmodel.hh
navierstokesncmodel.hh
volumevariables.hh
vtkoutputfields.hh
zeroeqncmodel.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/compositional)
