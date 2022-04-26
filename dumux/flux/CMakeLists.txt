add_subdirectory(box)
add_subdirectory(ccmpfa)
add_subdirectory(cctpfa)
add_subdirectory(porenetwork)
add_subdirectory(shallowwater)
add_subdirectory(staggered)

file(GLOB DUMUX_DUMUX_FLUX_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_FLUX_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/flux)
