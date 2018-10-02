add_subdirectory(binarycoefficients)
add_subdirectory(chemistry)
add_subdirectory(components)
add_subdirectory(constraintsolvers)
add_subdirectory(eos)
add_subdirectory(fluidmatrixinteractions)
add_subdirectory(fluidstates)
add_subdirectory(fluidsystems)
add_subdirectory(solidstates)
add_subdirectory(solidsystems)
add_subdirectory(spatialparams)

install(FILES
constants.hh
idealgas.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/material)
