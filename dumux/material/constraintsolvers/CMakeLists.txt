add_subdirectory(2p2c)

install(FILES
compositionalflash.hh
compositionfromfugacities.hh
computefromreferencephase.hh
immiscibleflash.hh
misciblemultiphasecomposition.hh
ncpflash.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material/constraintsolvers)
