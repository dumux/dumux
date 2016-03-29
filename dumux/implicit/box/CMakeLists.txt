
#install headers
install(FILES
intersectiontovertexbc.hh
assembler.hh
elementboundarytypes.hh
elementvolumevariables.hh
fvelementgeometry.hh
localresidual.hh
properties.hh
propertydefaults.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/box)
