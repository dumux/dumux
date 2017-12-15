add_subdirectory("omethod")

install(FILES
connectivitymap.hh
darcyslaw.hh
dualgridindexset.hh
elementfluxvariablescache.hh
elementvolumevariables.hh
fickslaw.hh
fluxvariablescachefiller.hh
fourierslaw.hh
fvelementgeometry.hh
fvgridgeometry.hh
generalconnectivitymap.hh
gridfluxvariablescache.hh
gridinteractionvolumeindexsets.hh
helper.hh
interactionvolumebase.hh
interactionvolumedatahandle.hh
interactionvolume.hh
methods.hh
properties.hh
subcontrolvolumeface.hh
tensorlambdafactory.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered/mpfa)
