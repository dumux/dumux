#install headers
install(FILES
fluxvariablescache.hh
indices.hh
model.hh
volumevariables.hh
fluxvariables.hh
localresidual.hh
problem.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/shallowwater/swe)
