add_subdirectory(sequential)

install(FILES
boxmaterialinterfaces.hh
formulation.hh
gridadaptindicator.hh
griddatatransfer.hh
incompressiblelocalresidual.hh
indices.hh
iofields.hh
model.hh
saturationreconstruction.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p)
