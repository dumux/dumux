add_subdirectory(solutions)

file(GLOB DUMUX_TEST_FREEFLOW_RANSNC_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_FREEFLOW_RANSNC_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/ransnc)
