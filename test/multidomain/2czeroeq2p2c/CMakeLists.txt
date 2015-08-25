add_input_file_links()

add_dumux_test(test_2czeroeq2p2c test_2czeroeq2p2c test_2czeroeq2p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_2czeroeq2p2c -ParameterFile ${CMAKE_CURRENT_SOURCE_DIR}/test_2czeroeq2p2c_reference.input"
                  --files ${CMAKE_SOURCE_DIR}/test/references/2czeroeq2p2c-ff-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2c-00028.vtu
                          ${CMAKE_SOURCE_DIR}/test/references/2czeroeq2p2c-pm-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/2p2c-00028.vtu)
