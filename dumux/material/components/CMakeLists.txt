add_subdirectory(iapws)

file(GLOB DUMUX_DUMUX_MATERIAL_COMPONENTS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_MATERIAL_COMPONENTS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/material/components)
