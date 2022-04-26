add_subdirectory(porenetwork)

file(GLOB DUMUX_DUMUX_MATERIAL_SPATIALPARAMS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_MATERIAL_SPATIALPARAMS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/material/spatialparams)
