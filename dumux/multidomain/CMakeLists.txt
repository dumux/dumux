add_subdirectory(boundary)
add_subdirectory(embedded)
add_subdirectory(facet)
add_subdirectory(io)

file(GLOB DUMUX_MULTIDOMAIN_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
