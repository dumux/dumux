add_subdirectory(cellcentered)
add_subdirectory(mimetic)
add_subdirectory(mpfa)

file(GLOB DUMUX_POROUSMEDIUMFLOW_2P_SEQUENTIAL_DIFFUSION_HEADERS *.hh *.inc)
install(FILES ${DUMUX_POROUSMEDIUMFLOW_2P_SEQUENTIAL_DIFFUSION_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/2p/sequential/diffusion)
