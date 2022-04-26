add_subdirectory(angeli)
add_subdirectory(channel)
add_subdirectory(donea)
add_subdirectory(kovasznay)
add_subdirectory(periodic)
add_subdirectory(sincos)

file(GLOB DUMUX_TEST_FREEFLOW_NAVIERSTOKES_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_FREEFLOW_NAVIERSTOKES_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/navierstokes)
