add_subdirectory("diffusion")
add_subdirectory("impes")
add_subdirectory("transport")

#install headers
install(FILES
indices.hh
properties.hh
celldataadaptive.hh
celldata.hh
fluxdata.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential)
