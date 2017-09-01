add_executable(test_2cstokes2p2c test_2cstokes2p2c test_2cstokes2p2c.cc)

# install sources
install(FILES
2cstokes2p2cproblem.hh
2cstokes2p2cspatialparams.hh
2p2csubproblem.hh
stokes2csubproblem.hh
test_2cstokes2p2c.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/2cstokes2p2c)