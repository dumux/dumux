add_dumux_test(test_zeroeq2c test_zeroeq2c test_zeroeq2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_zeroeq2c"
                  --files ${CMAKE_SOURCE_DIR}/test/references/zeroeq2c-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/zeroeq2c-00013.vtu)
