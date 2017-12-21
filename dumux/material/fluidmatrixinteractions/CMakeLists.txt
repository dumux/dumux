add_subdirectory("1p")
add_subdirectory("2p")
add_subdirectory("2pia")
add_subdirectory("3p")
add_subdirectory("mp")

#install headers
install(FILES
diffusivityconstanttortuosity.hh
diffusivitymillingtonquirk.hh
permeabilitykozenycarman.hh
permeabilityrutqvisttsang.hh
porosityprecipitation.hh
porosityreactivebed.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidmatrixinteractions)
