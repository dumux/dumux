file(GLOB DUMUX_POROUSMEDIUMFLOW_3PWATEROIL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_3PWATEROIL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/3pwateroil)
