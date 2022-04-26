add_subdirectory(staggered)

file(GLOB DUMUX_DUMUX_FREEFLOW_COMPOSITIONAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_FREEFLOW_COMPOSITIONAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/freeflow/compositional)
