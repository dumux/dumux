add_subdirectory("cellcentered")
add_subdirectory("mimetic")

#install headers
install(FILES
gridadapt.hh
gridadaptinitializationindicatordefault.hh
gridadaptinitializationindicator.hh
gridadaptproperties.hh
impet.hh
impetproblem.hh
impetproperties.hh
onemodelproblem.hh
pressureproperties.hh
properties.hh
transportproperties.hh
variableclassadaptive.hh
variableclass.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/sequential)
