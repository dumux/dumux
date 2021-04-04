add_subdirectory(cellcentered)
add_subdirectory(mimetic)

file(GLOB DUMUX_POROUSMEDIUMFLOW_SEQUENTIAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_SEQUENTIAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/sequential)
