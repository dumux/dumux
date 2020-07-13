add_subdirectory(1p)
add_subdirectory(1pnc)
add_subdirectory(2p)
add_subdirectory(common)

install(FILES
properties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetworkflow)
