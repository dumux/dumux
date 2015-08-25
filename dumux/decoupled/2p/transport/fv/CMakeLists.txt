
#install headers
install(FILES
capillarydiffusion.hh
convectivepart.hh
diffusivepart.hh
evalcflfluxcoats.hh
evalcflfluxdefault.hh
evalcflflux.hh
fvsaturation2p.hh
fvtransportproperties2p.hh
gravitypart.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/2p/transport/fv)
