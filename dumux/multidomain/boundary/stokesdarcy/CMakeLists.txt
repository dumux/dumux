file(GLOB DUMUX_MULTIDOMAIN_BOUNDARY_STOKESDARCY_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MULTIDOMAIN_BOUNDARY_STOKESDARCY_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/boundary/stokesdarcy)
