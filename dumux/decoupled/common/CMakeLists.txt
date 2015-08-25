add_subdirectory("fv")
add_subdirectory("mimetic")

#install headers
install(FILES
decoupledproperties.hh
gridadapt.hh
gridadaptinitializationindicatordefault.hh
gridadaptinitializationindicator.hh
gridadaptproperties.hh
impet.hh
impetproblem.hh
impetproperties.hh
onemodelproblem.hh
pressureproperties.hh
transportproperties.hh
variableclassadaptive.hh
variableclass.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/common)
