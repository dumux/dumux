add_subdirectory(convergencetest)

file(GLOB DUMUX_TEST_MULTIDOMAIN_BOUNDARY_STOKESDARCY_1P_1P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_MULTIDOMAIN_BOUNDARY_STOKESDARCY_1P_1P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/boundary/stokesdarcy/1p_1p)
