file(GLOB DUMUX_POROUSMEDIUMFLOW_CO2_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_CO2_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/co2)
