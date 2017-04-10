
#install headers
install(FILES
2pfluxvariables.hh
2plocalresidual.hh
2pmodel.hh
2pvolumevariables.hh
assembler.hh
elementvolumevariables.hh
fluxvariables.hh
localoperator.hh
model.hh
newtoncontroller.hh
propertydefaults.hh
simplecoupledproblem.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/geomechanics/el2p/decoupled)
