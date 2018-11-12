add_subdirectory(staggered)

install(FILES
fluxvariables.hh
indices.hh
kepsilonncmodel.hh
komegancmodel.hh
localresidual.hh
lowrekepsilonncmodel.hh
navierstokesncmodel.hh
oneeqncmodel.hh
volumevariables.hh
iofields.hh
zeroeqncmodel.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/compositional)
