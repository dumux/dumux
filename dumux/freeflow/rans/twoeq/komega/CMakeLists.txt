add_subdirectory(staggered)

file(GLOB DUMUX_DUMUX_FREEFLOW_RANS_TWOEQ_KOMEGA_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_FREEFLOW_RANS_TWOEQ_KOMEGA_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/freeflow/rans/twoeq/komega)
