add_subdirectory(oneeq)
add_subdirectory(twoeq)
add_subdirectory(zeroeq)

file(GLOB DUMUX_FREEFLOW_RANS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_RANS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans)
