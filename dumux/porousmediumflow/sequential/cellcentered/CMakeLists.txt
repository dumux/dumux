add_subdirectory(mpfa)

install(FILES
pressure.hh
transport.hh
velocity.hh
velocitydefault.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/porousmediumflow/sequential/cellcentered)
