add_subdirectory(1p)
add_subdirectory(1pia)
add_subdirectory(2p)
add_subdirectory(3p)
add_subdirectory(dispersiontensors)
add_subdirectory(frictionlaws)
add_subdirectory(mp)
add_subdirectory(porenetwork)

file(GLOB DUMUX_DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_HEADERS *.hh *.inc)
install(FILES ${DUMUX_DUMUX_MATERIAL_FLUIDMATRIXINTERACTIONS_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/dumux/material/fluidmatrixinteractions)
