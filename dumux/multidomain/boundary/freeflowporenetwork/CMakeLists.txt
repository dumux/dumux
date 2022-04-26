add_subdirectory(ffmassporenetwork)
add_subdirectory(ffmomentumporenetwork)

file(GLOB DUMUX_DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPORENETWORK_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOWPORENETWORK_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/multidomain/boundary/freeflowporenetwork)
