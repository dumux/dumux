
#install headers
install(FILES
stokesncfluxvariables.hh
stokesncindices.hh
stokesnclocalresidual.hh
stokesncmodel.hh
stokesncproperties.hh
stokesncpropertydefaults.hh
stokesncvolumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/stokesnc)
