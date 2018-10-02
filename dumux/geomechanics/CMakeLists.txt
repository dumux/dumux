add_subdirectory(elastic)
add_subdirectory(poroelastic)

install(FILES
fvproblem.hh
lameparams.hh
properties.hh
stressvariablescache.hh
velocityoutput.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/geomechanics)
