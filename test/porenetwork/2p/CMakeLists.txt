add_subdirectory(static)

file(GLOB DUMUX_TEST_PORENETWORK_2P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_PORENETWORK_2P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porenetwork/2p)
