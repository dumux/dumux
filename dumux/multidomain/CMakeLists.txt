add_subdirectory("2cnistokes2p2cni")
add_subdirectory("2cstokes2p2c")
add_subdirectory("common")
add_subdirectory("couplinglocalresiduals")

#install headers
install(FILES
assembler.hh
boxcouplinglocalresidual.hh
convergencewriter.hh
localoperator.hh
model.hh
newtoncontroller.hh
problem.hh
properties.hh
propertydefaults.hh
splitandmerge.hh
subdomainproperties.hh
subdomainpropertydefaults.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/multidomain)
