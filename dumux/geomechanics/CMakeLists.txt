add_subdirectory(elastic)
add_subdirectory(poroelastic)

file(GLOB DUMUX_DUMUX_GEOMECHANICS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_GEOMECHANICS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/geomechanics)
