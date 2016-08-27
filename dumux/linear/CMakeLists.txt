add_subdirectory(reordering)

#install headers
install(FILES
amgbackend.hh
amgparallelhelpers.hh
amgproperties.hh
linearsolverproperties.hh
pardisobackend.hh
seqsolverbackend.hh
vectorexchange.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/linear)
