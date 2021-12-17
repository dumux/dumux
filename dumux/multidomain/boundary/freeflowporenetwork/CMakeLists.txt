add_subdirectory(ffmassporenetwork)
add_subdirectory(ffmomentumporenetwork)

file(GLOB DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPORENETWORK_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPORENETWORK_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/boundary/freeflowporenetwork)
