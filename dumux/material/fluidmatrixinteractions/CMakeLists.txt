add_subdirectory("1p")
add_subdirectory("2p")
add_subdirectory("2pia")
add_subdirectory("3p")
add_subdirectory("mp")

#install headers
install(FILES
diffusivityconstanttau.hh
diffusivitymillingtonquirk.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidmatrixinteractions)
