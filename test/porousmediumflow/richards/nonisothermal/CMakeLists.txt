add_subdirectory(conduction)
add_subdirectory(convection)
add_subdirectory(evaporation)

file(GLOB DUMUX_TEST_POROUSMEDIUMFLOW_RICHARDS_NONISOTHERMAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_POROUSMEDIUMFLOW_RICHARDS_NONISOTHERMAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/richards/nonisothermal)
