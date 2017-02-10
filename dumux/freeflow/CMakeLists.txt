add_subdirectory("staggered")
add_subdirectory("stokes")
add_subdirectory("stokesnc")
add_subdirectory("stokesncni")
add_subdirectory("zeroeq")
add_subdirectory("zeroeqnc")
add_subdirectory("zeroeqncni")

#install headers
install(FILES
boundarylayermodel.hh
masstransfermodel.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow)
