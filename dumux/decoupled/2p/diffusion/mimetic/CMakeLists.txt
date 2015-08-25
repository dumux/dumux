
#install headers
install(FILES
croperator2padaptive.hh
croperator2p.hh
localstiffness.hh
mimetic2padaptive.hh
mimetic2p.hh
mimeticoperator2padaptive.hh
mimeticoperator2p.hh
mimeticpressure2padaptive.hh
mimeticpressure2p.hh
mimeticpressureproperties2padaptive.hh
mimeticpressureproperties2p.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/2p/diffusion/mimetic)
