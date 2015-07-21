add_dumux_test(test_el1p2c test_el1p2c test_el1p2c.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                 --script fuzzy
                 --files ${CMAKE_SOURCE_DIR}/test/references/el1p2c-reference.vtu
                         ${CMAKE_CURRENT_BINARY_DIR}/el1p2c-00015.vtu
                 --command "${CMAKE_CURRENT_BINARY_DIR}/test_el1p2c")
