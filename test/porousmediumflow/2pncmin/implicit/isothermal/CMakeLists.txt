dune_symlink_to_source_files(FILES params.input)

# isothermal tests
dumux_add_test(NAME test_2pncmin_dissolution_box
              LABELS porousmediumflow 2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=DissolutionBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncmin_dissolution_box-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_dissolution_box-00043.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_dissolution_box -ParameterFile params.input -Problem.Name test_2pncmin_dissolution_box")

dumux_add_test(NAME test_2pncmin_dissolution_box_restart
              LABELS porousmediumflow 2pncmin
              TARGET  test_2pncmin_dissolution_box
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS --script fuzzy
                       --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncmin_dissolution_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_dissolution_box_restart-00005.vtu
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_dissolution_box params.input -Problem.Name test_2pncmin_dissolution_box_restart -TimeLoop.DtInitial 50000 -Restart.Time 772823 -Restart.File test_2pncmin_dissolution_box-00038.vtu")

# the restart test has to run after the test that produces the corresponding vtu file
set_tests_properties(test_2pncmin_dissolution_box_restart PROPERTIES DEPENDS test_2pncmin_dissolution_box)

dumux_add_test(NAME test_2pncmin_dissolution_tpfa
              LABELS porousmediumflow 2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=DissolutionCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncmin_dissolution_tpfa-reference.vtu
                                    ${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_dissolution_tpfa-00042.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncmin_dissolution_tpfa -ParameterFile params.input -Problem.Name test_2pncmin_dissolution_tpfa")
