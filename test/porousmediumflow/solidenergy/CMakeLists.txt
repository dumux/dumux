add_input_file_links()

# non-isothermal tests
dumux_add_test(NAME test_solidenergy_tpfa
              LABELS porousmediumflow solidenergy
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_solidenergy_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_solidenergy_tpfa-00010.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_solidenergy_tpfa params.input -Problem.Name test_solidenergy_tpfa")
