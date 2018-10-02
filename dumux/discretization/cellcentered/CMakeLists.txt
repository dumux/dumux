add_subdirectory(mpfa)
add_subdirectory(tpfa)

install(FILES
connectivitymap.hh
elementboundarytypes.hh
elementsolution.hh
gridvolumevariables.hh
subcontrolvolume.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered)
