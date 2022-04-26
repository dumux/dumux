add_subdirectory(boundingboxtree)

file(GLOB DUMUX_TEST_GEOMETRY_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_GEOMETRY_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/geometry)
