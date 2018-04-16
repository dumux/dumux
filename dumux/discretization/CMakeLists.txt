add_subdirectory("box")
add_subdirectory("cellcentered")
add_subdirectory("staggered")

install(FILES
basefvgridgeometry.hh
darcyslaw.hh
evalgradients.hh
evalsolution.hh
fickslaw.hh
fluxstencil.hh
fluxvariablesbase.hh
fluxvariablescaching.hh
fourierslaw.hh
fvgridvariables.hh
fvproperties.hh
methods.hh
scvandscvfiterators.hh
scvoperator.hh
stationaryvelocityfield.hh
subcontrolvolumebase.hh
subcontrolvolumefacebase.hh
upwindscheme.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization)
