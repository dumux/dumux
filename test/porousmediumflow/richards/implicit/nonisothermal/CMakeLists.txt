add_subdirectory(conduction)
add_subdirectory(convection)
add_subdirectory(evaporation)

#install sources
install(FILES
spatialparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/implicit/richards/nonisothermal)
