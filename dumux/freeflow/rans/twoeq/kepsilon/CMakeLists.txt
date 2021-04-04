add_subdirectory(staggered)

file(GLOB DUMUX_FREEFLOW_RANS_TWOEQ_KEPSILON_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_RANS_TWOEQ_KEPSILON_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans/twoeq/kepsilon)
