add_subdirectory(static)

install(FILES
elementfluxvariablescache.hh
fluxvariablescache.hh
gridfluxvariablescache.hh
invasionstate.hh
iofields.hh
model.hh
newtonconsistencychecks.hh
newtonsolver.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porenetworkflow/2p)
