add_subdirectory(mpfa)
add_subdirectory(tpfa)

file(GLOB DUMUX_DUMUX_MULTIDOMAIN_FACET_CELLCENTERED_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_MULTIDOMAIN_FACET_CELLCENTERED_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/multidomain/facet/cellcentered)
