add_subdirectory("omethod")

install(FILES
computetransmissibility.hh
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
gridfluxvariablescache.hh
gridinteractionvolumeindexsets.hh
helper.hh
interactionvolumebase.hh
interactionvolumedatahandle.hh
localassembler.hh
localfacedata.hh
methods.hh
properties.hh
subcontrolvolumeface.hh
tensorlambdafactory.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/discretization/cellcentered/mpfa)
