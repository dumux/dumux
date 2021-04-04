add_subdirectory(sequential)

file(GLOB DUMUX_POROUSMEDIUMFLOW_2P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_2P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p)
