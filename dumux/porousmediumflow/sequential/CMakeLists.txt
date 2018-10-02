add_subdirectory(cellcentered)
add_subdirectory(mimetic)

install(FILES
gridadapt.hh
gridadaptinitializationindicator.hh
gridadaptinitializationindicatordefault.hh
gridadaptproperties.hh
impet.hh
impetproblem.hh
impetproperties.hh
onemodelproblem.hh
pressureproperties.hh
properties.hh
transportproperties.hh
variableclass.hh
variableclassadaptive.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/sequential)
