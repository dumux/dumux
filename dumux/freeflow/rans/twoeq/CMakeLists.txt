add_subdirectory(kepsilon)
add_subdirectory(komega)
add_subdirectory(lowrekepsilon)

file(GLOB DUMUX_FREEFLOW_RANS_TWOEQ_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_RANS_TWOEQ_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans/twoeq)
