add_subdirectory("diffusion")

#install headers
install(FILES
1pindices.hh
1pproperties.hh
celldata1p.hh
fluxdata1p.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/1p)
