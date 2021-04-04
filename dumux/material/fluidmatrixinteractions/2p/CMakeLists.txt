add_subdirectory(interfacialarea)
add_subdirectory(thermalconductivity)

file(GLOB DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_2P_HEADERS *.hh *.inc)
install(FILES ${DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_2P_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidmatrixinteractions/2p)
