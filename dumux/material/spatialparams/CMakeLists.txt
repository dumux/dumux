add_subdirectory(porenetwork)

file(GLOB DUMUX_MATERIAL_SPATIALPARAMS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MATERIAL_SPATIALPARAMS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/spatialparams)
