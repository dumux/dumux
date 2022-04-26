add_subdirectory(components)
add_subdirectory(fluidsystems)
add_subdirectory(solidsystems)

file(GLOB DUMUX_EXAMPLES_BIOMINERALIZATION_MATERIAL_HEADERS *.hh *.inc)
install(FILES ${DUMUX_EXAMPLES_BIOMINERALIZATION_MATERIAL_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/examples/biomineralization/material)
