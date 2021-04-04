add_subdirectory(box)
add_subdirectory(cellcentered)

file(GLOB DUMUX_MULTIDOMAIN_FACET_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_FACET_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/facet)
