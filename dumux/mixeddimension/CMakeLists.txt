add_subdirectory("embedded")
add_subdirectory("facet")
add_subdirectory("glue")

install(FILES
assembler.hh
integrationpointsource.hh
model.hh
newtoncontroller.hh
problem.hh
properties.hh
subproblemlocaljacobian.hh
subproblemproperties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/mixeddimension)
