add_subdirectory(staggered)

file(GLOB DUMUX_DUMUX_FREEFLOW_RANS_TWOEQ_SST_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_FREEFLOW_RANS_TWOEQ_SST_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/freeflow/rans/twoeq/sst)
