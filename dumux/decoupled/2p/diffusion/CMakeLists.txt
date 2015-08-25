add_subdirectory("fv")
add_subdirectory("fvmpfa")
add_subdirectory("mimetic")

#install headers
install(FILES
diffusionproblem2p.hh
diffusionproperties2p.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/2p/diffusion)
