
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
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/stokes)
