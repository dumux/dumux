
#install headers
install(FILES
boxassembler.hh
boxelementboundarytypes.hh
boxelementvolumevariables.hh
boxfvelementgeometry.hh
boxlocalresidual.hh
boxproperties.hh
boxpropertydefaults.hh
intersectiontovertexbc.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/box)
