
#install headers
install(FILES
assembler.hh
localjacobian.hh
localresidual.hh
model.hh
problem.hh
properties.hh
propertydefaults.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/fem)
