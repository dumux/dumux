add_subdirectory("twoeq")
add_subdirectory("zeroeq")

#install headers
install(FILES
model.hh
problem.hh
volumevariables.hh
vtkoutputfields.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/rans)
