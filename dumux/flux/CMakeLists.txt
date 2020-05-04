add_subdirectory(box)
add_subdirectory(ccmpfa)
add_subdirectory(cctpfa)
add_subdirectory(shallowwater)
add_subdirectory(staggered)

install(FILES
darcyslaw.hh
effectivestresslaw.hh
fickiandiffusioncoefficients.hh
fickslaw.hh
fluxvariablesbase.hh
fluxvariablescaching.hh
forchheimerslaw.hh
fourierslaw.hh
fourierslawnonequilibrium.hh
hookeslaw.hh
maxwellstefandiffusioncoefficients.hh
maxwellstefanslaw.hh
referencesystemformulation.hh
shallowwaterflux.hh
stationaryvelocityfield.hh
traits.hh
upwindscheme.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/flux)
