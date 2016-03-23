
#install headers
install(FILES
fluxvariables.hh
indices.hh
localjacobian.hh
localresidual.hh
model.hh
newtoncontroller.hh
problem.hh
properties.hh
propertydefaults.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/stokes)
