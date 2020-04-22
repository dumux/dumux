dune_symlink_to_source_files(FILES params.input)

# isothermal tests
dumux_add_test(NAME test_2pncminni_salinization_box
              LABELS porousmediumflow 2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=SalinizationBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncminni_salinization_box-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_box-00042.vtu
                            --zeroThreshold {"velocity_liq \(m/s\)_0":1e-18,"velocity_gas \(m/s\)_0":1e-14}
                               --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_box -ParameterFile params.input -Problem.Name test_2pncminni_salinization_box")

dumux_add_test(NAME test_2pncminni_salinization_tpfa
              LABELS porousmediumflow 2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=SalinizationCCTpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncminni_salinization_tpfa-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_tpfa-00042.vtu
                            --zeroThreshold {"velocity_liq \(m/s\)_0":1e-18,"velocity_gas \(m/s\)_0":1e-14}
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization_tpfa -ParameterFile params.input -Problem.Name test_2pncminni_salinization_tpfa")
