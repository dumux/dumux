
#install headers
install(FILES
ccassembler.hh
ccelementboundarytypes.hh
ccelementvolumevariables.hh
ccfvelementgeometry.hh
cclocalresidual.hh
ccproperties.hh
ccpropertydefaults.hh
assembler.hh
elementboundarytypes.hh
elementvolumevariables.hh
fvelementgeometry.hh
localresidual.hh
properties.hh
propertydefaults.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/cellcentered)
