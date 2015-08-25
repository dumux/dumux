add_subdirectory("fv")

#install headers
install(FILES
diffusionproblem1p.hh
diffusionproperties1p.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/1p/diffusion)
