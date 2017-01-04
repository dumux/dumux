add_subdirectory(adaptive)
add_subdirectory(box)
add_subdirectory(fem)
add_subdirectory(cellcentered)
add_subdirectory(cornerpoint)

install(FILES
assembler.hh
localjacobian.hh
localresidual.hh
model.hh
problem.hh
properties.hh
propertydefaults.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit)
