file(GLOB DUMUX_POROUSMEDIUMFLOW_2P2C_SEQUENTIAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_2P2C_SEQUENTIAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p2c/sequential)
