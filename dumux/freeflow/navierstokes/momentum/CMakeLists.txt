add_subdirectory(cvfe)
add_subdirectory(diamond)
add_subdirectory(pq1bubble)

file(GLOB DUMUX_FREEFLOW_NAVIERSTOKES_MOMENTUM_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_NAVIERSTOKES_MOMENTUM_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/navierstokes/momentum)
