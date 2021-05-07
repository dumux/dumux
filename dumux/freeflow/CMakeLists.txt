add_subdirectory(compositional)
add_subdirectory(navierstokes)
add_subdirectory(nonisothermal)
add_subdirectory(rans)
add_subdirectory(shallowwater)
add_subdirectory(advectiondiffusion)

file(GLOB DUMUX_FREEFLOW_HEADERS *.hh *.inc)
install(FILES ${DUMUX_FREEFLOW_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow)
