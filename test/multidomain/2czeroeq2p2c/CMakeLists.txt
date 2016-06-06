add_input_file_links()

add_dumux_test(test_2czeroeq2p2c test_2czeroeq2p2c test_2czeroeq2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2czeroeq2p2c -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2czeroeq2p2c_reference.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2czeroeq2p2c-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2c-00029.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2czeroeq2p2c-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2c-00029.vtu)

#install sources
install(FILES
2czeroeq2p2cproblem.hh
2czeroeq2p2cspatialparameters.hh
2p2csubproblem.hh
test_2czeroeq2p2c.cc
zeroeq2csubproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/multidomain/2czeroeq2p2c)
