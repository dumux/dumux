
#install headers
install(FILES
stokesncfluxvariables.hh
stokesncindices.hh
stokesnclocalresidual.hh
stokesncmodel.hh
stokesncproperties.hh
stokesncpropertydefaults.hh
stokesncvolumevariables.hh
fluxvariables.hh
indices.hh
localresidual.hh
model.hh
properties.hh
propertydefaults.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/stokesnc)
