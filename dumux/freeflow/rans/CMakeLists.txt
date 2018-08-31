add_subdirectory(oneeq)
add_subdirectory(twoeq)
add_subdirectory(zeroeq)

install(FILES
model.hh
problem.hh
volumevariables.hh
iofields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans)
