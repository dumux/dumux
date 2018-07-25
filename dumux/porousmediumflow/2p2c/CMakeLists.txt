add_subdirectory(sequential)

install(FILES
indices.hh
model.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p2c)
