add_subdirectory("box")
add_subdirectory("cellcentered")

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
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization)
