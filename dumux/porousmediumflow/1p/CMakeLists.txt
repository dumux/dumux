add_subdirectory(sequential)

install(FILES
incompressiblelocalresidual.hh
indices.hh
iofields.hh
model.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/1p)
