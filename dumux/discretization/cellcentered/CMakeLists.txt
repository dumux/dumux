add_subdirectory("tpfa")
add_subdirectory("mpfa")

install(FILES
connectivitymap.hh
elementboundarytypes.hh
elementsolution.hh
gridvolumevariables.hh
subcontrolvolume.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered)
