add_subdirectory(mpfa)
add_subdirectory(tpfa)
add_subdirectory(godunov)

install(FILES
connectivitymap.hh
elementboundarytypes.hh
elementsolution.hh
gridvolumevariables.hh
subcontrolvolume.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered)
