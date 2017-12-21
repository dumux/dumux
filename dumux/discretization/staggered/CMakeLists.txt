add_subdirectory("freeflow")

install(FILES
connectivitymap.hh
elementfacevariables.hh
elementfluxvariables.hh
elementvolumevariables.hh
facesolution.hh
fvelementgeometry.hh
fvgridgeometry.hh
globalfacevariables.hh
gridfluxvariables.hh
gridvariables.hh
gridvolumevariables.hh
properties.hh
subcontrolvolumeface.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/staggered)
