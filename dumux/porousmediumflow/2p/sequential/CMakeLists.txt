add_subdirectory(diffusion)
add_subdirectory(impes)
add_subdirectory(transport)

install(FILES
celldata.hh
celldataadaptive.hh
fluxdata.hh
indices.hh
properties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential)
