add_subdirectory("diffusion")
add_subdirectory("impes")
add_subdirectory("transport")

#install headers
install(FILES
2pindices.hh
2pproperties.hh
celldata2padaptive.hh
celldata2p.hh
fluxdata2p.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/2p)
