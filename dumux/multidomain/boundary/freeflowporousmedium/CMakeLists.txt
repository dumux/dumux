add_subdirectory(ffmasspm)
add_subdirectory(ffmomentumpm)

file(GLOB DUMUX_DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPOROUSMEDIUM_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPOROUSMEDIUM_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/multidomain/boundary/freeflowporousmedium)
