add_subdirectory(iapws)

file(GLOB DUMUX_MATERIAL_COMPONENTS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MATERIAL_COMPONENTS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/components)
