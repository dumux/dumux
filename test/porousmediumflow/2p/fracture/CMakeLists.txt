add_subdirectory(dt)
add_subdirectory(grids)

file(GLOB DUMUX_TEST_POROUSMEDIUMFLOW_2P_FRACTURE_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_POROUSMEDIUMFLOW_2P_FRACTURE_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/porousmediumflow/2p/fracture)
