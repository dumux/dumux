add_subdirectory(energy)
add_subdirectory(mass)
add_subdirectory(momentum)
add_subdirectory(staggered)

file(GLOB DUMUX_DUMUX_FREEFLOW_NAVIERSTOKES_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_FREEFLOW_NAVIERSTOKES_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/freeflow/navierstokes)
