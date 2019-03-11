add_subdirectory(box)
add_subdirectory(ccmpfa)
add_subdirectory(cctpfa)
add_subdirectory(shallowwater)
add_subdirectory(staggered)

install(FILES
darcyslaw.hh
effectivestresslaw.hh
fickslaw.hh
fluxvariablesbase.hh
fluxvariablescaching.hh
forchheimerslaw.hh
fourierslaw.hh
fourierslawnonequilibrium.hh
hookeslaw.hh
maxwellstefanslaw.hh
shallowwaterflux.hh
stationaryvelocityfield.hh
upwindscheme.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/flux)
