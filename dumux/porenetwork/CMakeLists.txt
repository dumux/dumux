add_subdirectory(1p)
add_subdirectory(1pnc)
add_subdirectory(2p)
add_subdirectory(2pnc)
add_subdirectory(common)
add_subdirectory(util)

file(GLOB DUMUX_PORENETWORK_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PORENETWORK_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetwork)
