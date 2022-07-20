add_subdirectory(box)
add_subdirectory(cellcentered)
add_subdirectory(facecentered)
add_subdirectory(fem)
add_subdirectory(porenetwork)
add_subdirectory(projection)
add_subdirectory(staggered)
add_subdirectory(cvfe)

file(GLOB DUMUX_DISCRETIZATION_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DISCRETIZATION_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization)
