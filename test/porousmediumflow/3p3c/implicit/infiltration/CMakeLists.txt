add_input_file_links(FILES params.input)

# isothermal tests
dumux_add_test(NAME test_3p3c_infiltration_box
              LABELS porousmediumflow 3p3c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePThreeCBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3c_infiltration_box-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_3p3c_infiltration_box-00005.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3c_infiltration_box params.input -Problem.Name test_3p3c_infiltration_box")

dumux_add_test(NAME test_3p3c_infiltration_tpfa
              LABELS porousmediumflow 3p3c
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=InfiltrationThreePThreeCCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_3p3c_infiltration_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_3p3c_infiltration_tpfa-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_3p3c_infiltration_tpfa params.input -Problem.Name test_3p3c_infiltration_tpfa")
