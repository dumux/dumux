add_subdirectory(localrules)
add_subdirectory(regularizedlocalrules)
install(FILES
baselocalrules.hh
localrules.hh
regularizedlocalrules.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/fluidmatrixinteractions/porenetwork/pore/2p)
