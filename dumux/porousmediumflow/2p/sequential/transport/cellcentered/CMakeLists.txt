
#install headers
install(FILES
capillarydiffusion.hh
convectivepart.hh
diffusivepart.hh
evalcflfluxcoats.hh
evalcflfluxdefault.hh
evalcflflux.hh
saturation.hh
transportproperties.hh
gravitypart.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential/transport/cellcentered)
