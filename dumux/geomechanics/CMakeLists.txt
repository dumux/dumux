add_subdirectory("elastic")

install(FILES
fvproblem.hh
lameparams.hh
properties.hh
stressvariablescache.hh
velocityoutput.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/geomechanics)
