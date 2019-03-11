#install headers
install(FILES
boundaryfluxes.hh
fluxvariables.hh
fluxvariablescache.hh
indices.hh
localresidual.hh
model.hh
properties.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflox/shallowwater/)
