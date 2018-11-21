add_subdirectory(staggered)

install(FILES
fluxvariables.hh
iofields.hh
kepsilonncmodel.hh
komegancmodel.hh
localresidual.hh
lowrekepsilonncmodel.hh
navierstokesncmodel.hh
oneeqncmodel.hh
volumevariables.hh
zeroeqncmodel.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/compositional)
