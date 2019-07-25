add_subdirectory(elastic)
add_subdirectory(poroelastic)

install(FILES
lameparams.hh
problem.hh
properties.hh
stressvariablescache.hh
velocityoutput.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/geomechanics)
