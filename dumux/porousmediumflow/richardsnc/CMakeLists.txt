file(GLOB DUMUX_POROUSMEDIUMFLOW_RICHARDSNC_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_RICHARDSNC_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/richardsnc)
