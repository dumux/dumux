add_subdirectory(1p)
add_subdirectory(1pia)
add_subdirectory(2p)
add_subdirectory(2pia)
add_subdirectory(3p)
add_subdirectory(frictionlaws)
add_subdirectory(mp)

install(FILES
diffusivityconstanttortuosity.hh
diffusivitymillingtonquirk.hh
fluidmatrixinteraction.hh
permeabilitykozenycarman.hh
porositydeformation.hh
porosityprecipitation.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidmatrixinteractions)
