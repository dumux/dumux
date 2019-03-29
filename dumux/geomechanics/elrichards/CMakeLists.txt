
#install headers
install(FILES
amgbackend.hh
assembler.hh
basemodel.hh
elementvolumevariables.hh
fluxvariables.hh
indices.hh
localjacobian.hh
localoperator.hh
localresidual.hh
model.hh
newtoncontroller.hh
properties.hh
propertydefaults.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/geomechanics/elrichards)
