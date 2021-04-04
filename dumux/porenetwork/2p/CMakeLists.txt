add_subdirectory(static)

file(GLOB DUMUX_PORENETWORK_2P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PORENETWORK_2P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetwork/2p)
