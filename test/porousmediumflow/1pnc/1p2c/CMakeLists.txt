add_subdirectory(isothermal)
add_subdirectory(nonisothermal)

file(GLOB DUMUX_TEST_POROUSMEDIUMFLOW_1PNC_1P2C_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_POROUSMEDIUMFLOW_1PNC_1P2C_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/1pnc/1p2c)
