add_subdirectory(elastic)
add_subdirectory(poroelastic)
add_subdirectory(stressstate)

file(GLOB DUMUX_GEOMECHANICS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_GEOMECHANICS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/geomechanics)
