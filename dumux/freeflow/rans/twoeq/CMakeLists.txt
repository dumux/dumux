add_subdirectory(kepsilon)
add_subdirectory(komega)
add_subdirectory(lowrekepsilon)
add_subdirectory(sst)

file(GLOB DUMUX_DUMUX_FREEFLOW_RANS_TWOEQ_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_FREEFLOW_RANS_TWOEQ_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/freeflow/rans/twoeq)
