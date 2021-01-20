add_input_file_links()

dumux_add_test(NAME test_2p2c_mpnc_comparison_box
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPTwoCComparisonBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2c_comparison_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_mpnc_comparison_box-00009.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_mpnc_comparison_box params.input -Problem.Name test_2p2c_mpnc_comparison_box")

dumux_add_test(NAME test_2p2c_mpnc_comparison_tpfa
              LABELS porousmediumflow 2p2c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TwoPTwoCComparisonCC
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_mpnc_2p2c_comparison_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_mpnc_comparison_tpfa-00009.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2p2c_mpnc_comparison_tpfa params.input -Problem.Name test_2p2c_mpnc_comparison_tpfa")
