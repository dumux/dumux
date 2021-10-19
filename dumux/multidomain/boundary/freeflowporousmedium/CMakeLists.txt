add_subdirectory(ffmasspm)
add_subdirectory(ffmomentumpm)

file(GLOB DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPOROUSMEDIUM_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPOROUSMEDIUM_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/boundary/freeflowporousmedium)
