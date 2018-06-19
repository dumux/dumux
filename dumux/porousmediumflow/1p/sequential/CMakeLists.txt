add_subdirectory(diffusion)

install(FILES
celldata.hh
fluxdata.hh
indices.hh
properties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/1p/sequential)
