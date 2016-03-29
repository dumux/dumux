
#install headers
install(FILES
localoperator.hh
newtoncontroller.hh
properties.hh
propertydefaults.hh
2p2ccouplinglocalresidual.hh
stokesnccouplinglocalresidual.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain/2cstokes2p2c)
