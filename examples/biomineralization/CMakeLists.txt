dune_symlink_to_source_files(FILES injection_checkpoints.dat injection_type.dat params.input)

# compile MICP simplified chemistry column setup
dumux_add_test(NAME example_biomineralization
               LABELS porousmediumflow example
               SOURCES main.cc
               COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
               CMD_ARGS --script fuzzy
                        --files ${CMAKE_SOURCE_DIR}/test/references/example_biomineralization-reference.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/example_biomineralization-00018.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/example_biomineralization params.input -Problem.Name  example_biomineralization")
