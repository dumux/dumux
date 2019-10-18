add_input_file_links("params.input")

dumux_add_test(NAME test_1p_pointsources_timedependent_tpfa
              LABELS porousmediumflow 1p
              SOURCES main.cc
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              COMPILE_DEFINITIONS TYPETAG=OnePSingularityTimeDependentCCTpfa
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_1p_pointsources_timedependent_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timedependent_tpfa-00003.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_1p_pointsources_timedependent_tpfa params.input -TimeLoop.TEnd 4
                       -Problem.Name test_1p_pointsources_timedependent_tpfa")
