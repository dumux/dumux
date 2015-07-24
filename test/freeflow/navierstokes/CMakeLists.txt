add_input_file_links()

add_dumux_test(test_navierstokes test_navierstokes test_navierstokes.cc
               python ${CMAKE_SOURCE_DIR}/bin/runtest.py
                  --script fuzzy
                  --command "${CMAKE_CURRENT_BINARY_DIR}/test_navierstokes"
                  --files ${CMAKE_SOURCE_DIR}/test/references/navierstokes-reference.vtu
                          ${CMAKE_CURRENT_BINARY_DIR}/navierstokes-00008.vtu)
