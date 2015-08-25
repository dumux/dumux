add_subdirectory("mpfa")

#install headers
install(FILES
fvpressure.hh
fvtransport.hh
fvvelocitydefault.hh
fvvelocity.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/decoupled/common/fv)
