add_subdirectory("mpfa")

#install headers
install(FILES
pressure.hh
transport.hh
velocitydefault.hh
velocity.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/sequential/cellcentered)
