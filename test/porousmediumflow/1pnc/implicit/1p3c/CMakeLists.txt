add_input_file_links()

# 1pnc tests
dune_add_test(NAME test_1pnc_maxwellstefan_tpfa
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MaxwellStefanOnePThreeCTestCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pnc_maxwellstefan_tpfa-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_maxwellstefan_tpfa-00032.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_maxwellstefan_tpfa params.input -Problem.Name test_1pnc_maxwellstefan_tpfa")

dune_add_test(NAME test_1pnc_maxwellstefan_box
              LABELS porousmediumflow 1p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MaxwellStefanOnePThreeCTestBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_1pnc_maxwellstefan_box-reference.vtu
                                 ${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_maxwellstefan_box-00091.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_1pnc_maxwellstefan_box params.input -Problem.Name test_1pnc_maxwellstefan_box")
