add_subdirectory("diffusion")

#install headers
install(FILES
indices.hh
properties.hh
celldata.hh
fluxdata.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/1p/sequential)
