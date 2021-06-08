add_subdirectory(grids)

file(GLOB DUMUX_TEST_IO_GRIDMANAGER_HEADERS *.hh *.inc)
install(FILES ${DUMUX_TEST_IO_GRIDMANAGER_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/io/gridmanager)
