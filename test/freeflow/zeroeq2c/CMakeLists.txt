add_input_file_links()

add_dumux_test(test_zeroeq2c test_zeroeq2c test_zeroeq2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq2c"
                  --files ${CMAKE_SOURCE_DIR}/test/references/zeroeq2c-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2c-00013.vtu)

#install sources
install(FILES
test_zeroeq2c.cc
zeroeq2ctestproblem.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/test/freeflow/zeroeq2c)
