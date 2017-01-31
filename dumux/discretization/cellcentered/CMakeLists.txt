add_subdirectory("tpfa")

install(FILES
globalvolumevariables.hh
subcontrolvolume.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered)
