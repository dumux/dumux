dune_symlink_to_source_files(FILES params.input)

# isothermal tests
dune_add_test(NAME test_2pncminni_salinization
              LABELS porousmediumflow 2pncmin
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=SalinizationBox
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS      --script fuzzy
                            --files ${CMAKE_SOURCE_DIR}/test/references/test_2pncminni_salinization-reference.vtu
                               ${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization-00034.vtu
                            --command "${CMAKE_CURRENT_BINARY_DIR}/test_2pncminni_salinization -ParameterFile params.input -Problem.Name test_2pncminni_salinization")
