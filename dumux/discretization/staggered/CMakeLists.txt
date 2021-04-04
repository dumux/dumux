add_subdirectory(freeflow)

file(GLOB DUMUX_DISCRETIZATION_STAGGERED_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DISCRETIZATION_STAGGERED_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/staggered)
