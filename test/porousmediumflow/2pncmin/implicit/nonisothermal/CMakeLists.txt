dune_symlink_to_source_files(FILES params.input)

# isothermal tests
dune_add_test(NAME test_2pncminni_salinization_box
              LABELS porousmediumflow 2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=SalinizationBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncminni_salinization_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_box-00035.vtu
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_box -ParameterFile params.input -Problem.Name test_2pncminni_salinization_box")

dune_add_test(NAME test_2pncminni_salinization_tpfa
              LABELS porousmediumflow 2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=SalinizationCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncminni_salinization_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_tpfa-00035.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_tpfa -ParameterFile params.input -Problem.Name test_2pncminni_salinization_tpfa")
