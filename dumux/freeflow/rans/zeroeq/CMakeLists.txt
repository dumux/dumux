file(GLOB DUMUX_FREEFLOW_RANS_ZEROEQ_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_RANS_ZEROEQ_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans/zeroeq)
