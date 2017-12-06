#install headers
install(FILES
amgbackend.hh
amgparallelhelpers.hh
amgtraits.hh
linearsolveracceptsmultitypematrix.hh
linearsolverproperties.hh
pardisobackend.hh
scotchbackend.hh
seqsolverbackend.hh
vectorexchange.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/linear)
