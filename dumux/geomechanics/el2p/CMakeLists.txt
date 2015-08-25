
#install headers
install(FILES
el2pamgbackend.hh
el2passembler.hh
el2pbasemodel.hh
el2pelementvolumevariables.hh
el2pfluxvariables.hh
el2pindices.hh
el2plocaljacobian.hh
el2plocaloperator.hh
el2plocalresidual.hh
el2pmodel.hh
el2pnewtoncontroller.hh
el2pproperties.hh
el2ppropertydefaults.hh
el2pvolumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/geomechanics/el2p)
