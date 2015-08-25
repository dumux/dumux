
#install headers
install(FILES
fvpressure2padaptive.hh
fvpressure2p.hh
fvpressureproperties2padaptive.hh
fvpressureproperties2p.hh
fvpressurevelocity2p.hh
fvpressurevelocityproperties2p.hh
fvvelocity2padaptive.hh
fvvelocity2p.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/2p/diffusion/fv)
