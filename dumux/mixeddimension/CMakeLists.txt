add_subdirectory("embedded")
add_subdirectory("glue")

install(FILES
assembler.hh
bulklocaljacobian.hh
lowdimlocaljacobian.hh
newtoncontroller.hh
problem.hh
model.hh
properties.hh
subproblemproperties.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/mixeddimension)
