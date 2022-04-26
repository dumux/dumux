add_subdirectory(format)
add_subdirectory(grid)
add_subdirectory(vtk)
add_subdirectory(xml)

file(GLOB DUMUX_DUMUX_IO_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_IO_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/io)
