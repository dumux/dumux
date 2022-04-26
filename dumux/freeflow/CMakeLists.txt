add_subdirectory(compositional)
add_subdirectory(navierstokes)
add_subdirectory(nonisothermal)
add_subdirectory(rans)
add_subdirectory(shallowwater)

file(GLOB DUMUX_DUMUX_FREEFLOW_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_FREEFLOW_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/freeflow)
