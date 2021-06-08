add_subdirectory(nonisothermal)

file(GLOB DUMUX_TEST_PORENETWORK_1P_NONCREEPINGFLOW_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_PORENETWORK_1P_NONCREEPINGFLOW_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porenetwork/1p/noncreepingflow)
