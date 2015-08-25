
#install headers
install(FILES
ccassembler.hh
ccelementboundarytypes.hh
ccelementvolumevariables.hh
ccfvelementgeometry.hh
cclocalresidual.hh
ccproperties.hh
ccpropertydefaults.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/cellcentered)
