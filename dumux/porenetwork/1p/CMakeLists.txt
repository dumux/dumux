file(GLOB DUMUX_PORENETWORK_1P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_PORENETWORK_1P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetwork/1p)
