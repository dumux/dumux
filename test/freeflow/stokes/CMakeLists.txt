add_input_file_links()

add_dumux_test(test_stokes test_stokes test_stokes.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_stokes"
                  --files ${CMAKE_SOURCE_DIR}/test/references/stokes-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/stokes-00014.vtu)
