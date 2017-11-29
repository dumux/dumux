add_input_file_links()

dune_add_test(NAME test_1pnc_box
              SOURCES test_1pnc.cc
              COMPILE_DEFINITIONS TYPETAG=OnePNCBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/1ptestcc-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/1p_tpfa-00010.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_box")
