add_input_file_links()

add_dumux_test(test_2cstokes2p2c test_2cstokes2p2c test_2cstokes2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2cstokes2p2c -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2cstokes2p2c_reference.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2cstokes2p2c-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes2c-00007.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2cstokes2p2c-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2c-00007.vtu)

#install sources
install(FILES
2cstokes2p2cproblem.hh
2cstokes2p2cspatialparams.hh
2p2csubproblem.hh
stokes2csubproblem.hh
test_2cstokes2p2c.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/2cstokes2p2c)
