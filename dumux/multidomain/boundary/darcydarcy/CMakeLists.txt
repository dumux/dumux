file(GLOB DUMUX_MULTIDOMAIN_BOUNDARY_DARCYDARCY_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_BOUNDARY_DARCYDARCY_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/boundary/darcydarcy)
