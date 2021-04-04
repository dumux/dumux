add_subdirectory(diffusion)
add_subdirectory(impes)
add_subdirectory(transport)

file(GLOB DUMUX_POROUSMEDIUMFLOW_2P_SEQUENTIAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_2P_SEQUENTIAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential)
