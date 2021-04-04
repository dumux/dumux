add_subdirectory(porenetwork)

file(GLOB DUMUX_IO_GRID_HEADERS *.hh *.inc)
install(FILES ${DUMUX_IO_GRID_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/io/grid)
