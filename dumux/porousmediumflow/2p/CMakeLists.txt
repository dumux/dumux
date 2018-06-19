add_subdirectory(sequential)

install(FILES
boxmaterialinterfaceparams.hh
formulation.hh
gridadaptindicator.hh
griddatatransfer.hh
incompressiblelocalresidual.hh
indices.hh
model.hh
saturationreconstruction.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p)
