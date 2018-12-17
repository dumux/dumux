add_input_file_links()

dune_add_test(NAME test_mpnc_2p2c_comparison_box
              SOURCES main.cc
              LABELS mpnc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2c_comparison_box-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_box-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_box params.input -Problem.Name test_mpnc_2p2c_comparison_box")

dune_add_test(NAME test_mpnc_2p2c_comparison_tpfa
              SOURCES main.cc
              LABELS mpnc
              COMPILE_DEFINITIONS TYPETAG=MPNCComparisonCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2c_comparison_tpfa-reference.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_tpfa-00009.vtu
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_mpnc_2p2c_comparison_tpfa params.input -Problem.Name test_mpnc_2p2c_comparison_tpfa")
