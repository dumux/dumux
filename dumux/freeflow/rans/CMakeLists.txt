add_subdirectory(oneeq)
add_subdirectory(twoeq)
add_subdirectory(zeroeq)

install(FILES
iofields.hh
model.hh
problem.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans)
