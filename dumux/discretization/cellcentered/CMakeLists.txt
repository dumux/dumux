add_subdirectory(mpfa)
add_subdirectory(tpfa)

file(GLOB DUMUX_DUMUX_DISCRETIZATION_CELLCENTERED_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_DISCRETIZATION_CELLCENTERED_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/discretization/cellcentered)
