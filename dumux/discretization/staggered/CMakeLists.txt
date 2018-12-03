add_subdirectory(freeflow)

install(FILES
elementfacevariables.hh
elementfluxvariablescache.hh
elementsolution.hh
facesolution.hh
fvelementgeometry.hh
fvgridgeometry.hh
gridfacevariables.hh
gridfluxvariablescache.hh
gridvariables.hh
subcontrolvolumeface.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/staggered)
