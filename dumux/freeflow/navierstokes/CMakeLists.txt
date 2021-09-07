add_subdirectory(mass)
add_subdirectory(momentum)
add_subdirectory(staggered)

file(GLOB DUMUX_FREEFLOW_NAVIERSTOKES_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_NAVIERSTOKES_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/navierstokes)
