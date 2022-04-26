add_subdirectory(box)
add_subdirectory(cellcentered)

file(GLOB DUMUX_DUMUX_MULTIDOMAIN_FACET_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_MULTIDOMAIN_FACET_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/multidomain/facet)
