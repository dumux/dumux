add_subdirectory(sequential)

file(GLOB DUMUX_POROUSMEDIUMFLOW_2P2C_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_2P2C_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p2c)
