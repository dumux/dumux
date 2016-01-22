
#install headers
install(FILES
stokesfluxvariables.hh
stokesindices.hh
stokeslocaljacobian.hh
stokeslocalresidual.hh
stokesmodel.hh
stokesnewtoncontroller.hh
stokesproblem.hh
stokesproperties.hh
stokespropertydefaults.hh
stokesvolumevariables.hh
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
