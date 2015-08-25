add_subdirectory("fv")

#install headers
install(FILES
transportproblem2p.hh
transportproperties2p.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/2p/transport)
